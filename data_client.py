from snowflake.sqlalchemy import URL
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy_mixins import AllFeaturesMixin


SF_USER = "RADAR_DEV_MRNA_LAB"
SF_PASSWORD = "HwSiHtGv5mVkG7GMFiRf"

engine = create_engine(
    URL(
        user=SF_USER,
        password=SF_PASSWORD,
        account="sanofi-emea_rnd",
        database="RADAR_DEV",
        schema="PSA_RESULTS_MRNA",
        warehouse="RADAR_DEV_WH_LAB",
    )
)

Base = declarative_base()
session = scoped_session(sessionmaker(bind=engine, autocommit=True))


class CustomBase(Base, AllFeaturesMixin):
    __abstract__ = True


CustomBase.set_session(session)






