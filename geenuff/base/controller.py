from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from . import orm

class Controller(object):
    """Base class for objects the manipulate a GeenuFF database."""

    def __init__(self, database_path, err_path='/dev/null'):
        self.database_path = database_path
        self.err_path = err_path
        self.engine, self.session = Controller.mk_session(self.database_path)

    @staticmethod
    def mk_session(database_path):
        engine = create_engine(Controller.full_db_path(database_path), echo=False)
        orm.Base.metadata.create_all(engine)
        session = sessionmaker(bind=engine)()
        return engine, session

    @staticmethod
    def full_db_path(path):
        sqlite_prefix = 'sqlite:///'
        if not path.startswith(sqlite_prefix):
            return sqlite_prefix + path
        return path

