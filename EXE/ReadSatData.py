#python -m py_compile ReadSatData_EMSA_MEDESS.py
def ReadSatData_EMSA_MEDESS(path_image,N_OS) :

    import sys
    import time
    from xml.dom.minidom import parse

#	path_image=sys.argv[1]
#	N_OS=sys.argv[2]

    TOT_index=int(float(N_OS))
    dom1=parse(path_image)

    time_in=dom1.getElementsByTagName("csn:timeStamp")[0].toxml()

    Coords=dom1.getElementsByTagName("csn:center")[0].getElementsByTagName("gml:Point")[0].getElementsByTagName("gml:pos")[0].toxml()

    split1=time_in.split('>')[1].split('<')[0].split('-')
    year=split1[0]
    month=split1[1]
    split2=split1[2].split('T')
    day=split2[0]
    split3=split2[1].split(':')
    hour=split3[0]
    minutes=split3[1]

    coord1=Coords.split('>')[1].split('<')[0].split(' ')
    Lat_dec=str(float(coord1[1]))
    Lon_dec=str(float(coord1[0]))
    Lat=coord1[1].split('.')
    Lat_degrees=Lat[0]
    Lat_min=(float(Lat_dec)-float(Lat_degrees))*60
    Lon=coord1[0].split('.')
    Lon_degrees=Lon[0]
    Lon_min=(float(Lon_dec)-float(Lon_degrees))*60

    document=year+"/"+month+"/"+day+" "+hour+":"+minutes+"   Date of satellite image"+"\n"

    SUM_points=0

    for index in range(TOT_index):
        Points=dom1.getElementsByTagName("csn:geometry")[index].getElementsByTagName("gml:posList")[0].toxml()
        points1=Points.split('>')[1].split('<')[0].split(' ')
        n_points0=len(points1)
        n_points=str(n_points0/2)
        SUM_points=SUM_points+n_points0/2

    document=document+str(SUM_points)+"        Number of data points"+"\n"
    document=document+"  lat     lon"+"\n"

    for index in range(TOT_index):
        Points=dom1.getElementsByTagName("csn:geometry")[index].getElementsByTagName("gml:posList")[0].toxml()
        points1=Points.split('>')[1].split('<')[0].split(' ')
        n_points0=len(points1)
        n_points=str(n_points0/2)

        for i in range(0, n_points0-1,2):
            points_split1=points1[i].split(',')
            points_lon=str(float(points_split1[0]))
            points_split2=points1[i+1].split(',')
            points_lat=str(float(points_split2[0]))
            document=document+points_lat+" "+points_lon+"\n"

    name="initial" 

    file1=file(name+".txt",'w')
    file1.writelines(document)
    file1.close()

    return (Lat_degrees,Lat_min,Lon_degrees,Lon_min)

