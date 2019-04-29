from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from . import orm

class Controller(object):
    """Base class for objects the manipulate a GeenuFF database."""

    def __init__(self, database_path, err_path='/dev/null'):
        self.database_path = Controller.full_db_path(database_path)
        self.err_path = err_path
        self.engine = None
        self.session = None
        self._mk_session()  # initializes session, engine, and insertion_queues

    def _mk_session(self):
        self.engine = create_engine(self.database_path, echo=False)
        orm.Base.metadata.create_all(self.engine)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()

    @staticmethod
    def full_db_path(path):
        sqlite_prefix = 'sqlite:///'
        if not path.startswith(sqlite_prefix):
            return sqlite_prefix + path
        return path

