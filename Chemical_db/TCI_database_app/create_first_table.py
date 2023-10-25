from sqlalchemy import *
from config import host, port, database, user, password

conn_str = f"postgresql://{user}:{password}@{host}/{database}"
engine = create_engine(conn_str)
connection = engine.connect()
metadata = MetaData()
first_tb = Table('first_table', metadata,
   Column("name", String(255), nullable=False),
   Column("CAS", String(255), nullable=False),
   Column("code", String(255), nullable=False),
   Column("product_number", String(255), primary_key=True),
   Column("cas_rn", String(255), nullable=False),
   Column("reaxys_registry_number", Integer, nullable=False),
   Column("pubchem_substance_id", Integer, nullable=False),
   Column("sdbs_aist_spectral_db", Integer, nullable=False),
   Column("merck_index_14", Integer, nullable=False),
   Column("mdl_number", String(255), nullable=False)
)
metadata.create_all(engine)
query = insert(first_tb).values(name='Abietic Acid',
                                CAS="514-10-3",
                                code='A0001',
                                product_number= 'a0001',
                                cas_rn= '514-10-3',
                                reaxys_registry_number= 2221451.0,
                                pubchem_substance_id= 87561707.0,
                                sdbs_aist_spectral_db= 1471.0,
                                merck_index_14= 7.0,
                                mdl_number= 'mfcd03423567')
ResultProxy = connection.execute(query)
