import sys

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from geenuff.base.orm import Coordinate, Genome, Feature
from geenuff.base.helpers import full_db_path


class ExportController(object):
    def __init__(self, db_path_in, with_features_only=True):
        self.db_path_in = db_path_in
        self._mk_session()
        if with_features_only:
            self.coordinate_query = self._coords_with_feature_query()
        else:
            self.coordinate_query = self._all_coords_query()

    def _mk_session(self):
        self.engine = create_engine(full_db_path(self.db_path_in), echo=False)
        self.session = sessionmaker(bind=self.engine)()

    def _check_genome_names(self, *argv):
        for names in argv:
            if names:
                genome_ids = self.session.query(Genome.id).filter(Genome.species.in_(names)).all()
                if len(genome_ids) != len(names):
                    print('One or more of the given genome names can not be found in the database')
                    exit()

    def _coords_with_feature_query(self):
        return self.session.query(Feature.coordinate_id).distinct()

    def _all_coords_query(self):
        return self.session.query(Coordinate.id)

    def _get_coords_by_genome(self, genomes, exclude):
        coordinate_ids_of_interest = self.coordinate_query()
        if genomes:
            print('Selecting the following genomes: {}'.format(genomes), file=sys.stderr)
            all_coord_ids = (self.session.query(Coordinate.id)
                             .join(Genome, Genome.id == Coordinate.genome_id)
                             .filter(Genome.species.in_(genomes))
                             .filter(Coordinate.id.in_(coordinate_ids_of_interest))
                             .all())
        else:
            if exclude:
                print('Selecting all genomes from {} except: {}'.format(self.db_path_in, exclude),
                      file=sys.stderr)
                all_coord_ids = (self.session.query(Coordinate.id)
                                 .join(Genome, Genome.id == Coordinate.genome_id)
                                 .filter(Genome.species.notin_(exclude))
                                 .filter(Coordinate.id.in_(coordinate_ids_of_interest))
                                 .all())
            else:
                print('Selecting all genomes from {}'.format(self.db_path_in),
                      file=sys.stderr)
                all_coord_ids = coordinate_ids_of_interest.all()

        return all_coord_ids
