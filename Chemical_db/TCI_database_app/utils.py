import psycopg2
from rdkit import Chem
import pandas as pd


class get_connection:

    def __init__(self, *args, **kwargs):

        self.args = args
        self.kwargs = kwargs
        self.conn = None
        self.cur = None

    def __enter__(self):

        self.conn = psycopg2.connect(*self.args, **self.kwargs)
        self.cur = self.conn.cursor()
        return self.conn, self.cur

    def __exit__(self, exc_type, exc_value, traceback):
        self.cur.close()
        self.conn.close()


def cancel_all(dbname, user, password, host, port='5432'):
    '''
    dbnameの有無は本質的には関係ないが、dbnameを与えない場合user名と同名のDBにアクセスしようとする。
    その際にそのDBが存在しない場合、アクセスに失敗するため、存在し、かつそのuserがアクセスできるDBを１つ適当に設定する必要がある。
    将来的には、なんとかする。
    '''

    with get_connection(dbname=dbname, user=user, password=password, host=host, port=port) as (conn, cur):
        '''pg_stat_activity view will have one row per server process,
        showing information related to the current activity of that process.
        https://www.postgresql.org/docs/current/monitoring-stats.html#MONITORING-PG-STAT-ACTIVITY-VIEW'''
        cur.execute('SELECT * FROM pg_stat_activity')
        df = pd.DataFrame(cur.fetchall(), columns=[col.name for col in cur.description])

    for pid in df[df['usename']==user]['pid']:
        with get_connection(dbname=dbname, user=user, password=password, host=host, port=port) as (conn, cur):
            '''Cancels a query. PG_CANCEL_BACKEND is functionally equivalent to
            the CANCEL command. You can cancel queries currently being run by your user.
            Superusers can cancel any query. pid = process ID of the query to be canceled.
            https://docs.aws.amazon.com/redshift/latest/dg/PG_CANCEL_BACKEND.html'''
            cur.execute(f'SELECT pg_cancel_backend({pid})')


def get_cansmi(smi):
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)
    except:
        return None
