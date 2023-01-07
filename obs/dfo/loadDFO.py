"""
These functions read bottle data from DFO.sqlite, which has CTD and bottle data
but there is no oxygen data for the CTD.

They also read in CTD data that does have oxygen from DFO_CTD.sqlite.
"""
import datetime as dt
import numpy as np
import pandas as pd
import gsw
import os
from sqlalchemy import create_engine, case
from sqlalchemy.orm import create_session 
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.sql import and_, or_, not_, func

def get_lims(datelims, latlims, lonlims):
    # if useable datelims, latlims, lonlims not provided, set broadly:
    if len(datelims)<2:
        datelims=(dt.datetime(1900,1,1),dt.datetime.now())
    ymd0 = (datelims[0].year, datelims[0].month, datelims[0].day)
    ymd1 = (datelims[1].year, datelims[1].month, datelims[1].day)
    if len(latlims)<2:
        latlims=(-91,91)
    if len(lonlims)<2:
        lonlims=(-361,361)
    return ymd0, ymd1, latlims, lonlims
    

def loadDFO_CTD(basedir='', dbname='DFO_CTD.sqlite',
    datelims=(),latlims=(),lonlims=(), xyt_only=False):
    """
    This loads DFO bottle data stored in an SQLite database.
        
    basedir is location of database
    dbname is database name
    datelims, if provided, loads only data between first and second datetime in tuple
        (second date not required)
    latlims, if provided, loads only data in range latlims[0]<=lat<latlims[1]
    lonlims, if provided, loads only data in range latlims[0]<=lat<latlims[1]
    """
    ymd0, ymd1, latlims, lonlims = get_lims(datelims, latlims, lonlims)
    
    # if db does not exist, exit
    dbpath=os.path.abspath(os.path.join(basedir,dbname))
    if not os.path.isfile(dbpath):
        raise Exception(f'ERROR: file {dbpath} not found')

    # load database structure and start session
    engine = create_engine('sqlite:///' + dbpath, echo = False)
    Base = automap_base()
    # reflect the tables in salish.sqlite:
    Base.prepare(engine, reflect=True)
    # mapped classes have been created
    # existing tables:
    StationTBL=Base.classes.StationTBL
    ObsTBL=Base.classes.ObsTBL
    CalcsTBL=Base.classes.CalcsTBL
    session = create_session(bind = engine, autocommit = False, autoflush = True)
    
    if xyt_only:
        qry=session.query(StationTBL.ID.label('Station'),
            StationTBL.StartYear.label('Year'),StationTBL.StartMonth.label('Month'),
            StationTBL.StartDay.label('Day'),StationTBL.StartHour.label('Hour'),
            StationTBL.Lat,StationTBL.Lon).\
            select_from(StationTBL).join(ObsTBL,ObsTBL.StationTBLID==StationTBL.ID).\
            join(CalcsTBL,CalcsTBL.ObsTBLID==ObsTBL.ID).\
            filter(and_(or_(StationTBL.StartYear>ymd0[0],
                and_(StationTBL.StartYear==ymd0[0], StationTBL.StartMonth>ymd0[1]),
                and_(StationTBL.StartYear==ymd0[0], StationTBL.StartMonth==ymd0[1], StationTBL.StartDay>=ymd0[2])),
                or_(StationTBL.StartYear<ymd1[0],
                and_(StationTBL.StartYear==ymd1[0],StationTBL.StartMonth<ymd1[1]),
                and_(StationTBL.StartYear==ymd1[0],StationTBL.StartMonth==ymd1[1], StationTBL.StartDay<ymd1[2])),
                StationTBL.Lat>=latlims[0],StationTBL.Lat<latlims[1],
                StationTBL.Lon>=lonlims[0],StationTBL.Lon<lonlims[1],
                StationTBL.Include==True,ObsTBL.Include==True,CalcsTBL.Include==True))
    else:
        SA=case([(CalcsTBL.Salinity_T0_C0_SA!=None, CalcsTBL.Salinity_T0_C0_SA)], else_=
            case([(CalcsTBL.Salinity_T1_C1_SA!=None, CalcsTBL.Salinity_T1_C1_SA)], else_=
            case([(CalcsTBL.Salinity_SA!=None, CalcsTBL.Salinity_SA)], else_= None)))
        CT=case([(CalcsTBL.Temperature_Primary_CT!=None, CalcsTBL.Temperature_Primary_CT)], else_=
            case([(CalcsTBL.Temperature_Secondary_CT!=None, CalcsTBL.Temperature_Secondary_CT)], else_=
            CalcsTBL.Temperature_CT))
        ZD=case([(ObsTBL.Depth!=None,ObsTBL.Depth)], else_= CalcsTBL.Z)
        FL=case([(ObsTBL.Fluorescence_URU_Seapoint!=None,ObsTBL.Fluorescence_URU_Seapoint)], else_=
            ObsTBL.Fluorescence_URU_Wetlabs)
        
        qry=session.query(StationTBL.ID.label('Station'),
            StationTBL.StartYear.label('Year'),StationTBL.StartMonth.label('Month'),
            StationTBL.StartDay.label('Day'),StationTBL.StartHour.label('Hour'),
            StationTBL.Lat,StationTBL.Lon,ZD.label('Z'),SA.label('SA'),CT.label('CT'),FL.label('Fluor'),
            ObsTBL.Oxygen_Dissolved_SBE.label('DO_mLL'),ObsTBL.Oxygen_Dissolved_SBE_1.label('DO_umolkg')).\
            select_from(StationTBL).join(ObsTBL,ObsTBL.StationTBLID==StationTBL.ID).\
            join(CalcsTBL,CalcsTBL.ObsTBLID==ObsTBL.ID).\
            filter(and_(or_(StationTBL.StartYear>ymd0[0],
                and_(StationTBL.StartYear==ymd0[0], StationTBL.StartMonth>ymd0[1]),
                and_(StationTBL.StartYear==ymd0[0], StationTBL.StartMonth==ymd0[1], StationTBL.StartDay>=ymd0[2])),
                or_(StationTBL.StartYear<ymd1[0],
                and_(StationTBL.StartYear==ymd1[0],StationTBL.StartMonth<ymd1[1]),
                and_(StationTBL.StartYear==ymd1[0],StationTBL.StartMonth==ymd1[1], StationTBL.StartDay<ymd1[2])),
                StationTBL.Lat>=latlims[0],StationTBL.Lat<latlims[1],
                StationTBL.Lon>=lonlims[0],StationTBL.Lon<lonlims[1],
                StationTBL.Include==True,ObsTBL.Include==True,CalcsTBL.Include==True))

    # write to pandas DataFrame
    df1=pd.read_sql_query(qry.statement, engine)
    if len(df1) == 0:
        return None
    else:
        if xyt_only:
            df1['dtUTC']=[dt.datetime(int(y),int(m),int(d))+dt.timedelta(hours=h) for y,m,d,h in
                zip(df1['Year'],df1['Month'],df1['Day'],df1['Hour'])]
            df1 = df1[['Station','Lon','Lat','dtUTC']]
        else:
            df1['dtUTC']=[dt.datetime(int(y),int(m),int(d))+dt.timedelta(hours=h)
                for y,m,d,h in zip(df1['Year'],df1['Month'],df1['Day'],df1['Hour'])]
            df1['Z'] = -df1['Z'] # fix sign of z to be positive up
            df1 = df1.drop(['Year','Month','Day','Hour'],axis=1)
        session.close()
        engine.dispose()
        return df1

def loadDFO_bottle(basedir='', dbname='DFO.sqlite',
        datelims=(),latlims=(),lonlims=(), xyt_only=False):
    """
    This loads DFO bottle data stored in an SQLite database.
        
    basedir is location of database
    dbname is database name
    datelims, if provided, loads only data between first and second datetime in tuple
        (second date not required)
    latlims, if provided, loads only data in range latlims[0]<=lat<latlims[1]
    lonlims, if provided, loads only data in range latlims[0]<=lat<latlims[1]
    """
    ymd0, ymd1, latlims, lonlims = get_lims(datelims, latlims, lonlims)

    # if db does not exist, exit
    dbpath=os.path.abspath(os.path.join(basedir,dbname))
    if not os.path.isfile(dbpath):
        raise Exception(f'ERROR: file {dbpath} not found')

    # load database structure and start session
    engine = create_engine('sqlite:///' + dbpath, echo = False)
    Base = automap_base()
    # reflect the tables in database:
    Base.prepare(engine, reflect=True)
    # mapped classes have been created:
    # existing tables -- assign the bottle sample related ones to variables:
    StationTBL=Base.classes.BOTStationTBL
    ObsTBL=Base.classes.BOTObsTBL
    CalcsTBL=Base.classes.BOTCalcsTBL
    session = create_session(bind = engine, autocommit = False, autoflush = True)

    if xyt_only:
        qry=session.query(StationTBL.ID.label('Station'),StationTBL.StartYear.label('Year'),
            StationTBL.StartMonth.label('Month'),StationTBL.StartDay.label('Day'),StationTBL.StartHour.label('Hour'),
            StationTBL.Lat,StationTBL.Lon).\
            select_from(StationTBL).join(ObsTBL,ObsTBL.StationTBLID==StationTBL.ID).\
            join(CalcsTBL,CalcsTBL.ObsID==ObsTBL.ID).\
            filter(and_(or_(StationTBL.StartYear>ymd0[0],
                and_(StationTBL.StartYear==ymd0[0], StationTBL.StartMonth>ymd0[1]),
                and_(StationTBL.StartYear==ymd0[0], StationTBL.StartMonth==ymd0[1], StationTBL.StartDay>=ymd0[2])),
                or_(StationTBL.StartYear<ymd1[0],
                and_(StationTBL.StartYear==ymd1[0],StationTBL.StartMonth<ymd1[1]),
                and_(StationTBL.StartYear==ymd1[0],StationTBL.StartMonth==ymd1[1], StationTBL.StartDay<ymd1[2])),
                StationTBL.Lat>=latlims[0],StationTBL.Lat<latlims[1],
                StationTBL.Lon>=lonlims[0],StationTBL.Lon<lonlims[1]))
    else:
        # In some cases, it useful to merge data from various columns 
        # (columns were created based on original field names to maintain consistency
        # with ASCII format data set)
        # TEOS-10 Absolute (technically Reference) Salinity (g/kg) previously calculated using
        # gsw package and stored in database
        SA=case([(CalcsTBL.Salinity_Bottle_SA!=None, CalcsTBL.Salinity_Bottle_SA)], else_=
            case([(CalcsTBL.Salinity_T0_C0_SA!=None, CalcsTBL.Salinity_T0_C0_SA)], else_=
            case([(CalcsTBL.Salinity_T1_C1_SA!=None, CalcsTBL.Salinity_T1_C1_SA)], else_=
            case([(CalcsTBL.Salinity_SA!=None, CalcsTBL.Salinity_SA)], else_=
            case([(CalcsTBL.Salinity__Unknown_SA!=None, CalcsTBL.Salinity__Unknown_SA)], 
            else_=CalcsTBL.Salinity__Pre1978_SA)))))
        CT=case([(CalcsTBL.Temperature_CT!=None, CalcsTBL.Temperature_CT)], else_=
            case([(CalcsTBL.Temperature_Primary_CT!=None, CalcsTBL.Temperature_Primary_CT)], else_=
            case([(CalcsTBL.Temperature_Secondary_CT!=None, CalcsTBL.Temperature_Secondary_CT)], 
            else_=CalcsTBL.Temperature_Reversing_CT)))
        # Define query to return desired data:
        # Labels will determine column headings in eventual pandas DataFrame
        #   - For columns from tables, format is tableVariable.ColumnName, and label is not required
        #   - For case statements defined above, label is always required
        # If data is to be returned from multiple tables, they must be joined
        # (see lines containing select_from and join below)
        # The 'Include' column was created for each table to flag rows with QC concerns; 
        # any rows with Include=False are not returned
        qry=session.query(StationTBL.ID.label('Station'),StationTBL.StartYear.label('Year'),
            StationTBL.StartMonth.label('Month'),StationTBL.StartDay.label('Day'),StationTBL.StartHour.label('Hour'),
            StationTBL.Lat,StationTBL.Lon,ObsTBL.Pressure,ObsTBL.Depth,
            ObsTBL.Chlorophyll_Extracted.label('Chl'),ObsTBL.Chlorophyll_Extracted_units.label('Chl_units'),
            ObsTBL.Nitrate_plus_Nitrite.label('N'),ObsTBL.Nitrate_plus_Nitrite_units.label('N_units'),
            ObsTBL.Ammonium.label('NH4'),ObsTBL.Ammonium_units.label('NH4_units'),
            ObsTBL.Silicate.label('SiO4'),ObsTBL.Silicate_units.label('SiO4_units'),
            SA.label('SA'),CT.label('CT'),
            ObsTBL.Oxygen_Dissolved.label('DO'),ObsTBL.Oxygen_Dissolved_units.label('DO_units')).\
            select_from(StationTBL).join(ObsTBL,ObsTBL.StationTBLID==StationTBL.ID).\
            join(CalcsTBL,CalcsTBL.ObsID==ObsTBL.ID).\
            filter(and_(or_(StationTBL.StartYear>ymd0[0],
                and_(StationTBL.StartYear==ymd0[0], StationTBL.StartMonth>ymd0[1]),
                and_(StationTBL.StartYear==ymd0[0], StationTBL.StartMonth==ymd0[1], StationTBL.StartDay>=ymd0[2])),
                or_(StationTBL.StartYear<ymd1[0],
                and_(StationTBL.StartYear==ymd1[0],StationTBL.StartMonth<ymd1[1]),
                and_(StationTBL.StartYear==ymd1[0],StationTBL.StartMonth==ymd1[1], StationTBL.StartDay<ymd1[2])),
                StationTBL.Lat>=latlims[0],StationTBL.Lat<latlims[1],
                StationTBL.Lon>=lonlims[0],StationTBL.Lon<lonlims[1]))

    # write to pandas DataFrame
    df1 = pd.read_sql(qry.statement, engine)
    # print(df1.columns)
    if (len(df1) == 0) or ('Depth' not in df1.columns):
        return None
    else:
        if xyt_only:
            df1['dtUTC']=[dt.datetime(int(y),int(m),int(d))+dt.timedelta(hours=h) for y,m,d,h in
                zip(df1['Year'],df1['Month'],df1['Day'],df1['Hour'])]
            df1 = df1[['Station','Lon','Lat','dtUTC']]
        else:
            # Z will be positive up
            df1.loc[df1.Depth.isnull(), 'Depth'] = np.nan
            if 'Pressure' in df1.columns:
                df1.loc[df1.Pressure.isnull(), 'Pressure'] = np.nan
                df1['Z']=np.where(df1.Depth >= 0,-df1['Depth'],gsw.z_from_p(p=df1['Pressure'].to_numpy(),lat=df1['Lat'].to_numpy()))
            else:
                df1['Z'] = -df1['Depth']
            df1['dtUTC']=[dt.datetime(int(y),int(m),int(d))+dt.timedelta(hours=h) for y,m,d,h in
                zip(df1['Year'],df1['Month'],df1['Day'],df1['Hour'])]
            if df1.Z.isna().all():
                df1 = None
            else:
                df1 = df1.drop(['Pressure','Depth','Year','Month','Day','Hour'],axis=1)
        session.close()
        engine.dispose()
        return df1
    

