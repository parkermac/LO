"""
These functions read bottle data from DFO.sqlite, whic has CTD and bottle data
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

def loadDFO_CTD(basedir='', dbname='DFO_CTD.sqlite', datelims=()):
    """
    load DFO CTD data stored in SQLite database (exclude most points outside Salish Sea)
    basedir is location of database
    dbname is database name
    datelims, if provided, loads only data between first and second datetime in tuple
    """
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
    SA=case([(CalcsTBL.Salinity_T0_C0_SA!=None, CalcsTBL.Salinity_T0_C0_SA)], else_=
             case([(CalcsTBL.Salinity_T1_C1_SA!=None, CalcsTBL.Salinity_T1_C1_SA)], else_=
             case([(CalcsTBL.Salinity_SA!=None, CalcsTBL.Salinity_SA)], else_= None)))
    CT=case([(CalcsTBL.Temperature_Primary_CT!=None, CalcsTBL.Temperature_Primary_CT)], else_=
             case([(CalcsTBL.Temperature_Secondary_CT!=None, CalcsTBL.Temperature_Secondary_CT)], else_=CalcsTBL.Temperature_CT))
    ZD=case([(ObsTBL.Depth!=None,ObsTBL.Depth)], else_= CalcsTBL.Z)
    FL=case([(ObsTBL.Fluorescence_URU_Seapoint!=None,ObsTBL.Fluorescence_URU_Seapoint)], else_= ObsTBL.Fluorescence_URU_Wetlabs)
    if len(datelims)<2:
        qry=session.query(StationTBL.StartYear.label('Year'),StationTBL.StartMonth.label('Month'),
                      StationTBL.StartDay.label('Day'),StationTBL.StartHour.label('Hour'),
                      StationTBL.Lat,StationTBL.Lon,ZD.label('Z'),SA.label('SA'),CT.label('CT'),FL.label('Fluor'),
                      ObsTBL.Oxygen_Dissolved_SBE.label('DO_mLL'),ObsTBL.Oxygen_Dissolved_SBE_1.label('DO_umolkg')).\
                select_from(StationTBL).join(ObsTBL,ObsTBL.StationTBLID==StationTBL.ID).\
                join(CalcsTBL,CalcsTBL.ObsTBLID==ObsTBL.ID).filter(and_(StationTBL.Lat>47-3/2.5*(StationTBL.Lon+123.5),
                    StationTBL.Lat<47-3/2.5*(StationTBL.Lon+121),
                    StationTBL.Include==True,ObsTBL.Include==True,CalcsTBL.Include==True))
    else:
        start_y=datelims[0].year
        start_m=datelims[0].month
        start_d=datelims[0].day
        end_y=datelims[1].year
        end_m=datelims[1].month
        end_d=datelims[1].day
        qry=session.query(StationTBL.ID.label('CTDStationTBLID'),StationTBL.StartYear.label('Year'),StationTBL.StartMonth.label('Month'),
                      StationTBL.StartDay.label('Day'),StationTBL.StartHour.label('Hour'),
                      StationTBL.Lat,StationTBL.Lon,ZD.label('Z'),SA.label('SA'),CT.label('CT'),FL.label('Fluor'),
                      ObsTBL.Oxygen_Dissolved_SBE.label('DO_mLL'),ObsTBL.Oxygen_Dissolved_SBE_1.label('DO_umolkg')).\
                select_from(StationTBL).join(ObsTBL,ObsTBL.StationTBLID==StationTBL.ID).\
                join(CalcsTBL,CalcsTBL.ObsTBLID==ObsTBL.ID).filter(and_(or_(StationTBL.StartYear>start_y,
                         and_(StationTBL.StartYear==start_y, StationTBL.StartMonth>start_m),
                         and_(StationTBL.StartYear==start_y, StationTBL.StartMonth==start_m, StationTBL.StartDay>=start_d)),
                        or_(StationTBL.StartYear<end_y,
                         and_(StationTBL.StartYear==end_y,StationTBL.StartMonth<end_m),
                         and_(StationTBL.StartYear==end_y,StationTBL.StartMonth==end_m, StationTBL.StartDay<end_d)),
                    StationTBL.Lat>47-3/2.5*(StationTBL.Lon+123.5),
                    StationTBL.Lat<47-3/2.5*(StationTBL.Lon+121),
                    StationTBL.Include==True,ObsTBL.Include==True,CalcsTBL.Include==True))
    df1=pd.read_sql_query(qry.statement, engine)
    df1['dtUTC']=[dt.datetime(int(y),int(m),int(d))+dt.timedelta(hours=h) for y,m,d,h in zip(df1['Year'],df1['Month'],df1['Day'],df1['Hour'])]
    df1['Z'] = -df1['Z'] # fix sign of z to be positive up
    session.close()
    engine.dispose()
    return df1

def loadDFO_bottle(basedir='.', dbname='DFO.sqlite',
        datelims=(),latlims=(),lonlims=(), xyt_only=False):
    """
    This loads bottle data.
        
    load DFO data stored in SQLite database
    basedir is location of database
    dbname is database name
    datelims, if provided, loads only data between first and second datetime in tuple (second
      date not included)
    latlims, if provided, loads only data in range latlims[0]<=lat<latlims[1]
    lonlims, if provided, loads only data in range latlims[0]<=lat<latlims[1]
    lonlims
    """
    # definitions
    # if useable datelims, latlims, lonlims not provided, set broadly:
    if len(datelims)<2:
        datelims=(dt.datetime(1900,1,1),dt.datetime.now())
    start_y=datelims[0].year
    start_m=datelims[0].month
    start_d=datelims[0].day
    end_y=datelims[1].year
    end_m=datelims[1].month
    end_d=datelims[1].day
    if len(latlims)<2:
        latlims=(-91,91)
    if len(lonlims)<2:
        lonlims=(-361,361)

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
            filter(and_(or_(StationTBL.StartYear>start_y,
                and_(StationTBL.StartYear==start_y, StationTBL.StartMonth>start_m),
                and_(StationTBL.StartYear==start_y, StationTBL.StartMonth==start_m, StationTBL.StartDay>=start_d)),
                or_(StationTBL.StartYear<end_y,
                and_(StationTBL.StartYear==end_y,StationTBL.StartMonth<end_m),
                and_(StationTBL.StartYear==end_y,StationTBL.StartMonth==end_m, StationTBL.StartDay<end_d)),
                StationTBL.Lat>=latlims[0],
                StationTBL.Lat<latlims[1],
                StationTBL.Lon>=lonlims[0],
                StationTBL.Lon<lonlims[1]))
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
        # If data is to be returned from multiple tables, they must be joined (see lines containing select_from and join below)
        # The 'Include' column was created for each table to flag rows with QC concerns; 
        # any rows with Include=False are not returned
        qry=session.query(StationTBL.ID.label('Station'),StationTBL.StartYear.label('Year'),
            StationTBL.StartMonth.label('Month'),StationTBL.StartDay.label('Day'),StationTBL.StartHour.label('Hour'),
            StationTBL.Lat,StationTBL.Lon,#ObsTBL.sourceFile,
            ObsTBL.Pressure,ObsTBL.Depth,
            ObsTBL.Chlorophyll_Extracted.label('Chl'),ObsTBL.Chlorophyll_Extracted_units.label('Chl_units'),
            ObsTBL.Nitrate_plus_Nitrite.label('N'),ObsTBL.Nitrate_plus_Nitrite_units.label('N_units'),
            ObsTBL.Silicate.label('Si'),ObsTBL.Silicate_units.label('Si_units'),
            SA.label('SA'),CT.label('CT'),
            ObsTBL.Oxygen_Dissolved.label('DO'),ObsTBL.Oxygen_Dissolved_units.label('DO_units')).\
            select_from(StationTBL).join(ObsTBL,ObsTBL.StationTBLID==StationTBL.ID).\
            join(CalcsTBL,CalcsTBL.ObsID==ObsTBL.ID).\
            filter(and_(or_(StationTBL.StartYear>start_y,
                and_(StationTBL.StartYear==start_y, StationTBL.StartMonth>start_m),
                and_(StationTBL.StartYear==start_y, StationTBL.StartMonth==start_m, StationTBL.StartDay>=start_d)),
                or_(StationTBL.StartYear<end_y,
                and_(StationTBL.StartYear==end_y,StationTBL.StartMonth<end_m),
                and_(StationTBL.StartYear==end_y,StationTBL.StartMonth==end_m, StationTBL.StartDay<end_d)),
                StationTBL.Lat>=latlims[0],
                StationTBL.Lat<latlims[1],
                StationTBL.Lon>=lonlims[0],
                StationTBL.Lon<lonlims[1]))

    # write to pandas DataFrame
    df1 = pd.read_sql(qry.statement, engine)
    if xyt_only:
        df1['dtUTC']=[dt.datetime(int(y),int(m),int(d))+dt.timedelta(hours=h) for y,m,d,h in zip(df1['Year'],df1['Month'],df1['Day'],df1['Hour'])]
        df1 = df1[['Station','Lon','Lat','dtUTC']]
    else:
        df1['Z']=np.where(df1['Depth']>=0,df1['Depth'],gsw.z_from_p(p=df1['Pressure'].to_numpy(),lat=df1['Lat'].to_numpy()))
        df1['dtUTC']=[dt.datetime(int(y),int(m),int(d))+dt.timedelta(hours=h) for y,m,d,h in zip(df1['Year'],df1['Month'],df1['Day'],df1['Hour'])]
        df1 = df1.drop(['Pressure','Depth','Year','Month','Day','Hour'],axis=1)
    session.close()
    engine.dispose()
    return df1
    
def loadDFO_CTD_orig(basedir='.', dbname='DFO.sqlite',
        datelims=(),latlims=(),lonlims=()):
    """
    This version works with the original database DFO.sqlite, which has both bottle and
    CTD data, but the CTD data dod not have oxygen. I keep it only for completeness.
        
    load DFO CTD data stored in SQLite database (exclude most points outside Salish Sea)
    basedir is location of database
    dbname is database name
    datelims, if provided, loads only data between first and second datetime in tuple (second
      date not included)
    latlims, if provided, loads only data in range latlims[0]<=lat<latlims[1]
    lonlims, if provided, loads only data in range latlims[0]<=lat<latlims[1]
    lonlims
    """
    # definitions
    # if useable datelims not provided, set broadly:
    if len(datelims)<2:
        datelims=(dt.datetime(1900,1,1),dt.datetime.now())
    start_y=datelims[0].year
    start_m=datelims[0].month
    start_d=datelims[0].day
    end_y=datelims[1].year
    end_m=datelims[1].month
    end_d=datelims[1].day
    if len(latlims)<2:
        latlims=(-91,91)
    if len(lonlims)<2:
        lonlims=(-361,361)

    # if db does not exist, exit
    dbpath=os.path.abspath(os.path.join(basedir,dbname))
    if not os.path.isfile(dbpath):
        raise Exception(f'ERROR: file {dbpath} not found')

    # load database structure and start session
    engine = create_engine('sqlite:///' + dbpath, echo = False)
    Base = automap_base()
    # reflect the tables in database:
    Base.prepare(engine, reflect=True)
    # mapped classes have been created for
    # existing tables -- assign the CTD-related ones to variables:
    StationTBL=Base.classes.CTDStationTBL
    ObsTBL=Base.classes.CTDObsTBL
    CalcsTBL=Base.classes.CTDCalcsTBL
    session = create_session(bind = engine, autocommit = False, autoflush = True)

    # In some cases, it useful to merge data from variaous columns 
    #  (columns were created based on original field names to maintain consistency with ASCII format data set)
    # TEOS-10 Absolute (technically Reference) Salinity (g/kg) previously calculated using gsw package and stored in database
    SA=case([(CalcsTBL.Salinity_T0_C0_SA!=None, CalcsTBL.Salinity_T0_C0_SA)], else_=
             case([(CalcsTBL.Salinity_T1_C1_SA!=None, CalcsTBL.Salinity_T1_C1_SA)], else_=
             case([(CalcsTBL.Salinity_SA!=None, CalcsTBL.Salinity_SA)], else_= None)))
    # TEOS-10 Conservative Temperature (deg C) previously calculated using gsw package and stored in database
    CT=case([(CalcsTBL.Temperature_Primary_CT!=None, CalcsTBL.Temperature_Primary_CT)], else_=
             case([(CalcsTBL.Temperature_Secondary_CT!=None, CalcsTBL.Temperature_Secondary_CT)], 
                   else_=CalcsTBL.Temperature_CT))
    # If a depth in meters was associated with the original data, return that; otherwise,
    #  return the stored depth calculated from pressure and latitude using gsw
    ZD=case([(ObsTBL.Depth!=None,ObsTBL.Depth)], else_= CalcsTBL.Z)
    # Fluorescence: return data from either type of instrument, choosing Seapoint where available 
    # (majority of measurements available are Seapoint)
    FL=case([(ObsTBL.Fluorescence_URU_Seapoint!=None,ObsTBL.Fluorescence_URU_Seapoint)], 
             else_= ObsTBL.Fluorescence_URU_Wetlabs)

    # Define query to return desired data:
    # Labels will determine column headings in eventual pandas DataFrame
    #   - For columns from tables, format is tableVariable.ColumnName, and label is not required
    #   - For case statements defined above, label is always required
    # If data is to be returned from multiple tables, they must be joined (see lines containing select_from and join below)
    # The 'Include' column was created for each table to flag rows with QC concerns; 
    #   any rows with Include=False are not returned
    qry=session.query(StationTBL.ID.label('CTDStationTBLID'),StationTBL.StartYear.label('Year'),StationTBL.StartMonth.label('Month'),
                  StationTBL.StartDay.label('Day'),StationTBL.StartHour.label('Hour'),
                  StationTBL.Lat,StationTBL.Lon,ZD.label('Z'),SA.label('SA'),CT.label('CT'),FL.label('Fluor')).\
            select_from(StationTBL).join(ObsTBL,ObsTBL.StationTBLID==StationTBL.ID).\
            join(CalcsTBL,CalcsTBL.ObsTBLID==ObsTBL.ID).\
            filter(and_(or_(StationTBL.StartYear>start_y,
                            and_(StationTBL.StartYear==start_y, StationTBL.StartMonth>start_m),
                            and_(StationTBL.StartYear==start_y, StationTBL.StartMonth==start_m, StationTBL.StartDay>=start_d)),
                        or_(StationTBL.StartYear<end_y,
                            and_(StationTBL.StartYear==end_y,StationTBL.StartMonth<end_m),
                            and_(StationTBL.StartYear==end_y,StationTBL.StartMonth==end_m, StationTBL.StartDay<end_d)),
                        StationTBL.Lat>=latlims[0],
                        StationTBL.Lat<latlims[1],
                        StationTBL.Lon>=lonlims[0],
                        StationTBL.Lon<lonlims[1],
                        StationTBL.Include==True,
                        ObsTBL.Include==True,
                        CalcsTBL.Include==True))
    df1 = pd.read_sql(qry.statement, engine)
    # df1=pd.DataFrame(qry.all())
    df1['dtUTC']=[dt.datetime(int(y),int(m),int(d))+dt.timedelta(hours=h) for y,m,d,h in zip(df1['Year'],df1['Month'],df1['Day'],df1['Hour'])]
    session.close()
    engine.dispose()
    return df1