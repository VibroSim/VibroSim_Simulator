import collections
import numpy as np


def read(filename):
    fh=open(filename,"r")

    metadata=collections.OrderedDict()
    fieldheaderdata = collections.OrderedDict()
    fieldrowdata = []
    rownum=0

    lines=[]
    line=fh.readline()
    while line != "":
        lines.append(line)
        line=fh.readline()
        pass

    independentvar = None # extra independent variable (different columns)
    
    # first read general metadata
    for linenum in range(len(lines)):
        line = lines[linenum]
        splitline =  line.split()
        #if splitline[0]=="%" and splitline[1].endswith(":"):
        #if line[0]=="%" and line.find(":") >= 0:
        if line[0:2]=="% " and linenum < len(lines)-1 and lines[linenum+1][0:2]=="% " and line.find(":") >= 0:
            # This is a general metadata line... we can tell because the following line still begins with '%'
            colonpos=line.find(":") # find first colon
            valuepos=colonpos+1
            while line[valuepos].isspace():
                valuepos+=1
                pass
            name=line[2:valuepos].strip()
            value=line[valuepos:].strip()
            metadata[splitline[1][:-1]]=value
            pass

        elif splitline[0]=="%":
            # Actual data header... this is a bunch of column headers
            fieldpos=1
            numfields=0
            while fieldpos < len(splitline):
                fieldname=splitline[fieldpos]
                units=None
                independentval=None
                pointcoords=[]
                
                #print("fieldpos=%d" % (fieldpos))
                fieldpos+=1

                
                while fieldpos < len(splitline)-1 and splitline[fieldpos]=='+':

                    # multiple terms being added together
                    fieldname = fieldname + " " + (" ".join(splitline[fieldpos:(fieldpos+2)]))
                    fieldpos += 2
                    pass
                
                    
                if fieldpos < len(splitline) and splitline[fieldpos].startswith('(') and splitline[fieldpos].endswith(')'):
                    # found bare units in parentheses after fieldname
                    units=splitline[fieldpos][1:-1]
                    fieldpos+=1
                    pass

                if fieldpos < len(splitline)-1 and splitline[fieldpos].startswith('(') and splitline[fieldpos].endswith('),') and splitline[fieldpos+1] == "Point:":
                    # found "(units), Point:(... )" after fieldname
                    
                    units=splitline[fieldpos][1:-2]
                    fieldpos+=2

                    # Now extract the point coordinates
                    startflag=True
                    endflag=False
                    while not endflag:
                        pointcoord=splitline[fieldpos]
                        fieldpos+=1
                        if startflag:
                            assert(pointcoord.startswith("("))
                            pointcoord=pointcoord[1:]
                            startflag=False
                            pass
                        
                        if pointcoord.find(')') >= 0:
                            # close parentheses means point coords are done
                            endflag=True
                            endpos=pointcoord.find(')')
                            # At least some versions of COMSOL leave no space
                            # before the next entry... so we break this
                            # splitline entry into two at endpos+1
                            nextlinepart = pointcoord[endpos+1:]
                            pointcoord=pointcoord[:endpos]

                            if len(nextlinepart) > 0:
                                splitline.insert(fieldpos,nextlinepart)
                                pass
                            pass
                        else:
                            # Coordinates other than the last should end with a comma that we need to strip
                            assert(pointcoord.endswith(","))
                            pointcoord=pointcoord[:-1]
                            pass
                        pointcoords.append(pointcoord)
                        pass
                    pass

                
                if fieldpos < len(splitline)-1 and splitline[fieldpos]=="@" and splitline[fieldpos+1].find('=') >= 0:
                    # got independent variable such as time
                    new_independentvar = splitline[fieldpos+1][:splitline[fieldpos+1].find('=')]
                    if independentvar is not None:
                        assert(independentvar == new_independentvar) # Can only handle one such independent variable
                        pass
                    else:
                        independentvar = new_independentvar
                        pass
                    fieldpos+=1
                    independentval = float(splitline[fieldpos][(len(independentvar)+1):])  # often time
                    fieldpos+=1
                    pass
                
                fullfieldname=fieldname
                if units is not None:
                    fullfieldname="%s (%s)" % (fieldname,units)
                    pass
                if len(pointcoords) > 0:
                    fullfieldname+=", Point: (" + ", ".join(pointcoords) + ")"
                    pass

                
                if independentval is None:
                    assert(fullfieldname not in fieldheaderdata) # Field can only be present once, unless at multiple times or (or other independent value) 
                    fieldheaderdata[fullfieldname]=numfields # numfields represents column # for following data
                    pass
                else:
                    if fullfieldname not in fieldheaderdata:
                        fieldheaderdata[fullfieldname]=[]
                        pass
                    
                    fieldheaderdata[fullfieldname].append((independentvar,independentval,numfields)) # numfields represents column # for following data 
                    pass
                
                numfields+=1

                
                pass
            pass
        else:
            # actual data
            assert(len(splitline)==numfields) # should be exactly one column per field
            fieldrowdata.append(collections.OrderedDict())
            
            for fieldname in list(fieldheaderdata.keys()):
                if isinstance(fieldheaderdata[fieldname],int):
                    # a column number
                    colnum=fieldheaderdata[fieldname]
                    try:
                        fieldrowdata[rownum][fieldname]=float(splitline[colnum])
                        pass
                    except ValueError:
                        fieldrowdata[rownum][fieldname]=complex(splitline[colnum].replace('i','j'))
                        pass
                    pass
                else:
                    # Iterate over list of (time,colnum)
                    independentvalarray=np.zeros(len(fieldheaderdata[fieldname]),dtype='d')
                    dataarray=np.zeros(len(fieldheaderdata[fieldname]),dtype=np.complex128)
                    index=0
                    for (independentvar_junk,independentval,colnum) in fieldheaderdata[fieldname]:
                        independentvalarray[index]=independentval
                        dataarray[index]=complex(splitline[colnum].replace('i','j'))
                        
                        index+=1
                        pass

                    if (dataarray.imag==0.0).all():
                        dataarray=dataarray.real
                        pass
                    
                    fieldrowdata[rownum][fieldname]=(independentvalarray,dataarray)
                    pass
                pass
            rownum+=1
            pass
        pass
    return (metadata,fieldheaderdata,fieldrowdata)
    
if __name__=="__main__":
    from function_as_script import scriptify

    read_script=scriptify(read)

    inputfile = "vibrodynamic_xducercontactprobe.txt"
    #inputfile="cantilever_model_probe_table_2019-06-15.txt"
    
    (metadata,fieldheaderdata,fieldrowdata)=read_script(inputfile)

    from matplotlib import pyplot as pl

    #pl.figure(1)
    #pl.clf()
    #pl.plot(fieldrowdata[0]["vibrodynamic_xducercontactprobe_displ (m)"][0],fieldro#wdata[0]["vibrodynamic_xducercontactprobe_displ (m)"][1])
    #pl.xlabel('Time (s)')

    col0 = np.array([fieldrowdatum[fieldheaderdata.keys()[0]] for fieldrowdatum in fieldrowdata])
    col1 = np.array([fieldrowdatum[fieldheaderdata.keys()[1]] for fieldrowdatum in fieldrowdata])
    pl.show()
