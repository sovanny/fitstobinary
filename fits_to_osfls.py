import math, numpy, struct, ctypes, time, shutil, collections
from os import walk
from enum import Enum
from astropy.io import fits
from astropy.time import Time

rootpath = "/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/"
writepath = "./"

def populatefiles(folder, fileending, daterange, type):
    for (dirpath, dirnames, filenames) in walk(rootpath + folder):
        for file in filenames:    
            if(file.endswith(fileending)):
                date = file[daterange[0]:daterange[1]]
                if(type == 'PFSS_IO'): fitsfilesdates[date] = []
                try: fitsfilesdates[date].append((folder + '/' + file, type))
                except: continue       
        break

# POPULATE THE FILE INDEX VARIABLE
fitsfilesdates = {} 
populatefiles("trace_pfss_intoout",'fits.gz', (0,12), "PFSS_IO")
populatefiles("trace_pfss_outtoin", 'fits.gz', (0,12), "PFSS_OI")
populatefiles("trace_scs_outtoin", 'fits.gz', (0,12), "SCS_OI")
populatefiles("WSA_OUT", 'fits', (4,16), "WSA_OUT")
populatefiles("WSA_VEL", 'fits', (4,16), "WSA_VEL")

datestodelete = []
for date in fitsfilesdates:
    if(len(fitsfilesdates[date]) != 5):
        datestodelete.append(date)
for date in datestodelete:
    del fitsfilesdates[date]
print("Deleted {} dates because they were incomplete sets".format(len(datestodelete)))
nDates = len(fitsfilesdates)
print("Total {} dates to convert.".format(nDates))

ordereddates = collections.OrderedDict(sorted(fitsfilesdates.items()))

# FUNCTION DEFINITIONS

def sph2cart(coord):
    return [
         coord[2] * math.sin(coord[1]) * math.cos(coord[0]), 
         coord[2] * math.sin(coord[1]) * math.sin(coord[0]),
         coord[2] * math.cos(coord[1])
    ]

def coord2datarow(coord):
    sph = sph2cart(coord[0:3])
    sph.append(coord[5])
    return [numpy.float32(x) for x in sph] 

class Model(Enum):
    Batsrus = 0
    Enlil = 1
    Pfss = 2
    Wsa = 3
    Invalid = 5

def obstimeToJ2000(ot):
    timefits = ot[0:4] + '-' + ot[5:7] + '-' + ot[8:10] + 'T'
    timefits += ot[11:13] + ot[14:17] + ot[18:21] + '.000'
    pathSafeTimeString = timefits.replace(':', '-')
    time = Time(timefits, format='fits')
    y2000 = Time(2000, format='jyear')
    jdaysdelta = time.jd - y2000.jd
    jsecs = jdaysdelta*60*60*24
    return [jsecs, pathSafeTimeString]

def toOsfls(filename, modelname, indices, typename, extra_indices, extra_file):
    fl_fits = fits.open(rootpath + filename)
    fl_data = fl_fits[0].data
    fl_fits.close()
    
    if(modelname != 'PFSS_IO'):
        vel_fits = fits.open(rootpath + extra_file)
        vel_data = vel_fits[0].data
        vel_data = vel_data[1] # the second layer has the wind speed
        vel_fits.close()
        indices.extend(extra_indices)
        indices = list(numpy.unique(indices))
        
    versionNumber = 0
    [triggerTime, pathSafeTimeString] = obstimeToJ2000(fl_fits[0].header['OBSTIME'])
    #print("TriggerTime: ",triggerTime)
    #print("pathSafeTimeString: ",pathSafeTimeString)
    fileName = pathSafeTimeString + '.osfls'

    model = Model.Wsa.value
    isMorphable = False

    nVert = 0
    lineStart = []
    lineCount = []
    vertexPositions = []

    extraQuantities_1 = []
    extraQuantities_2= []

    for i in indices:
        points = [coord2datarow(pt) for pt in fl_data[i] if pt[0] > -900] 
        if (len(points) < 2): continue
        lineStart.append(nVert)
        nVert += len(points)
        lineCount.append(len(points))
        [vertexPositions.extend(pt[0:3]) for pt in points] # extend to unfold elements
        
        if(modelname == 'PFSS_IO'):
            fieldstrength = [pt[3] for pt in points]
            extraQuantities_1.extend(fieldstrength)
            if(i in extra_indices): openness = 0
            elif(points[0][3] < 0): openness = 1 # open and negative
            else: openness = 2 # open and positive
            xtra = [openness]*len(points) # same value for all points
            extraQuantities_2.extend(xtra)
        else:
            in_sheet = 1 # 1 if not in sheet
            if(points[0][3] < 0): in_sheet = -1 # -1 if not in sheet and negative
            if(i in extra_indices): in_sheet = 2 # 2 if in sheet
            if(points[0][3] < 0): in_sheet = -2 # 2 if in sheet and negative
            xtra = [in_sheet]*len(points) # same value for all points
            extraQuantities_1.extend(xtra)
            velocity = vel_data[math.floor(i/180)][i%180] # speed at 21.5
            if(points[0][3] < 0): velocity = (-1)*velocity # negative polarity
            xtra = [velocity]*len(points) # same value for all points
            extraQuantities_2.extend(xtra)
    
    nLines = len(lineStart)

    nExtras = 2
    extraQuantities = extraQuantities_1
    extraQuantities.extend(extraQuantities_2)
   
    if(modelname == 'PFSS_IO'):
        extraQuantityNames = ['Field strength \0', 'Open/closed regions \0']
    else: 
        extraQuantityNames = ['Polarity Sheet \0', 'Solar wind speed and Polarity\ 0']
    nStringBytes = sum([len(s) for s in extraQuantityNames])
    allNamesInOne = ''
    for s in extraQuantityNames:
        allNamesInOne += s
    
    # Prepare data for writing to binary. Using Struct and pack
    typestr = '= i d i ? Q Q Q Q %sl %sL %sf %sf %ss' % (nLines, nLines, 3*nVert, nExtras*nVert, nStringBytes)
    struct_to_write = struct.Struct(typestr)
    #print('Uses           :', struct_to_write.size, 'bytes')
    values_to_write = (versionNumber, triggerTime, model, isMorphable, nLines, nVert, nExtras, nStringBytes)
    values_to_write += (*lineStart, *lineCount, *vertexPositions, *extraQuantities, allNamesInOne.encode('utf-8'))
    buffer = ctypes.create_string_buffer(struct_to_write.size)  
    struct_to_write.pack_into(buffer, 0, *values_to_write)
    fout = open(writepath + modelname + '/' + fileName, 'wb')
    fout.write(buffer)
    fout.close()
# end def toOsfls

def everynth(step):
    return (range(0,16200, step), 'step' + str(step))

def pickClosed(filename):
    fl_fits = fits.open(rootpath + filename)
    fl_data = fl_fits[0].data
    indices_to_save = []

    for i in range(16200):
        b_radii = [pt[5] for pt in fl_data[i] if pt[0] > -900] 
        if (len(b_radii) < 2): continue
        first_b_radius = b_radii[0]
        last_b_radius = b_radii[-1]
        if(first_b_radius*last_b_radius < 0): #if product is negative/opposite signs
            indices_to_save.append(i)
    return indices_to_save

def make_sparser(indices, step):
    new_indices = indices[::3]
    return [new_indices, 'closed_step' + str(step)]

def pick_sheet_lines(wsa_file):
    wsa_fits = fits.open(wsa_file)
    wsa_data = wsa_fits[0].data
    wsa_data = wsa_data[0] # the first layer has the coronal field information
    wsa_fits.close()

    sheet_indices = []
    for i in range(90):
        for j in range(180):
            if(i < 89 and numpy.sign(wsa_data[i,j]) != numpy.sign(wsa_data[i+1][j])):
                index = i*179 + j
                sheet_indices.append(index)
                index = (i+1)*179 + j
                sheet_indices.append(index)
            if(j < 179 and numpy.sign(wsa_data[i,j]) != numpy.sign(wsa_data[i][j+1])):
                index = i*179 + j
                sheet_indices.append(index)
                index = i*179 + (j+1)
                sheet_indices.append(index)

    return list(numpy.unique(sheet_indices))




def boundaryLines(filename_io, filename_oi, indices_closed):
    fl_io = (fits.open(rootpath + filename_io))[0].data
    fl_oi = (fits.open(rootpath + filename_oi))[0].data

    threshold = math.sqrt(2)
    # boundary_lines_io = [] do we want these?
    boundary_lines_oi = []
    io_last_phi = {}
    io_last_theta = {}
    io_first_phi = {}
    io_first_theta = {}

    for j in indices_closed:
        last = -1 # find last index
        for point in fl_io[j]:
            if (point[0] > -900):
                last += 1
            else: break
        # get the first and last coordinates for in-to-out (both on surface)
        io_last_phi[j] = math.degrees(fl_io[j][last][0])
        io_last_theta[j] = math.degrees(fl_io[j][last][1])
        io_first_phi[j] = math.degrees(fl_io[j][0][0])
        io_first_theta[j] = math.degrees(fl_io[j][0][1])

    for i in range(16200):
        last = -1
        for point in fl_oi[i]:
            if (point[0] > -900):
                last += 1
            else: break
        if (last < 2): continue
        # get the last coordinates for out-to-in (the ones on the surface)
        oi_last_phi = math.degrees(fl_oi[i][last][0])
        oi_last_theta = math.degrees(fl_oi[i][last][1])
        for j in indices_closed:
            # calculate distances and compare
            dist_last_phi = abs(io_last_phi[j] - oi_last_phi)
            dist_last_theta = abs(io_last_theta[j] - oi_last_theta)
            dist_first_phi = abs(io_first_phi[j] - oi_last_phi)
            dist_first_theta = abs(io_first_theta[j] - oi_last_theta)
            if(dist_last_phi < threshold and dist_last_theta < threshold):
                #boundary_lines_io.append(j)
                boundary_lines_oi.append(i)
                break
            if(dist_first_phi < threshold and dist_first_theta < threshold):
                #boundary_lines_io.append(j)
                boundary_lines_oi.append(i)
                break
    return [boundary_lines_oi, 'boundary']
# end def pickBoundaryLines

def fillOut(boundary_indices, step):
    output = []
    previous_index = 0
    for index in boundary_indices:
        subrange = list(range(previous_index, index, step))
        output.extend(subrange)
        if(output[-1] != index):    
            output.append(index)
        previous_index = index
    return [output, 'boundary_filled']

def mergeAndConvertToOsfls(filenamePFSS, filenameSCS):
    fl_fits_pfss = fits.open(rootpath + filenamePFSS)
    fl_data_pfss = fl_fits_pfss[0].data
    fl_fits_pfss.close()
    
    fl_fits_scs = fits.open(rootpath + filenameSCS)
    fl_data_scs = fl_fits_scs[0].data
    fl_fits_scs.close()
    
    indices = range(16200,16740)
  
    [triggerTime, pathSafeTimeString] = obstimeToJ2000(fl_fits_pfss[0].header['OBSTIME'])
    [triggerTime_s, pathSafeTimeString_s] = obstimeToJ2000(fl_fits_scs[0].header['OBSTIME'])
    if(triggerTime != triggerTime_s):
        print("Files are from different times. Canceling conversion. triggertimes: ")
        print(triggerTime)
        print(triggerTime_s)
        
    versionNumber = 0
    fileName = pathSafeTimeString + '.osfls'
    model = Model.Wsa.value
    isMorphable = False
    
    nVert = 0
    lineStart = []
    lineCount = []
    vertexPositions = []
    extraQuantities = []
    
    extraQuantities_polarity = []
    extraQuantities_level = []

    for i in indices:
        points_pfss = [coord2datarow(pt) for pt in fl_data_pfss[i] if pt[0] > -900] 
        points_scs = [coord2datarow(pt) for pt in fl_data_scs[i] if pt[0] > -900] 
        
        if (len(points_pfss) < 2 or len(points_scs) < 2): continue
            
        lineStart.append(nVert)
        
        combinedVertLen = (len(points_pfss) + len(points_scs))
        nVert += combinedVertLen
        lineCount.append(combinedVertLen)
        
        [vertexPositions.extend(pt[0:3]) for pt in points_scs] # add scs vertices FIRST
        [vertexPositions.extend(pt[0:3]) for pt in points_pfss] # add pfss vertices
  
        #polarity = 0
        #if(points_scs[0][3] >= 0):
        #    polarity = 1
        #xtra = [polarity]*combinedVertLen
        #extraQuantities_polarity.extend(xtra)

        fieldstrength = [pt[3] for pt in points_scs]
        extraQuantities_polarity.extend(fieldstrength)
        fieldstrength = [pt[3] for pt in points_pfss]
        extraQuantities_polarity.extend(fieldstrength)
        
        sunearthpoint = 1
        if (i in range(16380, 16560)):
            sunearthpoint = 0
        elif (i in range(16560, 16740)):
            sunearthpoint = 2
        xtra = [sunearthpoint]*combinedVertLen
        extraQuantities_level.extend(xtra)
       
    nLines = len(lineStart)

    nExtras = 2
    extraQuantities = extraQuantities_polarity
    extraQuantities.extend(extraQuantities_level)

    extraQuantityNames = ['Field strength \0', 'Sun-earth level\0']

    nStringBytes = sum([len(s) for s in extraQuantityNames])
    allNamesInOne = ''
    for s in extraQuantityNames:
        allNamesInOne += s
    
    # Prepare data for writing to binary. Using Struct and pack
    typestr = '= i d i ? Q Q Q Q %sl %sL %sf %sf %ss' % (nLines, nLines, 3*nVert, nExtras*nVert, nStringBytes)
    struct_to_write = struct.Struct(typestr)
    #print('Format string  :', struct_to_write.format)
    #print('Uses           :', struct_to_write.size, 'bytes')
    values_to_write = (versionNumber, triggerTime, model, isMorphable, nLines, nVert, nExtras, nStringBytes)
    values_to_write += (*lineStart, *lineCount, *vertexPositions, *extraQuantities, allNamesInOne.encode('utf-8'))
    
    buffer = ctypes.create_string_buffer(struct_to_write.size)    
    struct_to_write.pack_into(buffer, 0, *values_to_write)
    
    fout = open(writepath + 'SUN_EARTH/' + fileName, 'wb')
    fout.write(buffer)
    fout.close()

# ACTUAL PROGRAM
start_time = time.time()
count = 0
# relies heavily on having a sorted list of files for each timestep
# so that pfss_io goes first and picks out the lines
for timestamp in fitsfilesdates:
    pfss_io = fitsfilesdates[timestamp][0]
    pfss_oi = fitsfilesdates[timestamp][1]
    scs_oi = fitsfilesdates[timestamp][2]
    wsafile = fitsfilesdates[timestamp][3]
    velfile = fitsfilesdates[timestamp][4]
    indices_closed_lines = pickClosed(pfss_io[0])
    #indices = make_sparser(indices_closed_lines, 3) # change here for sparser
    indices = everynth(7)
    toOsfls(pfss_io[0], pfss_io[1], indices[0], indices[1], indices_closed_lines, 0)
    #print('Finished converting {} after {} seconds.'.format(pfss_io[1],time.time()-start_time))
    indices = boundaryLines(pfss_io[0], pfss_oi[0], indices_closed_lines)
    indices = fillOut(indices[0], 25) # change step here for sparser
    sheet_indices = pick_sheet_lines(rootpath + wsafile[0])
    toOsfls(pfss_oi[0], pfss_oi[1], indices[0], indices[1], sheet_indices, velfile[0])
    toOsfls(scs_oi[0], scs_oi[1], indices[0], indices[1], sheet_indices, velfile[0])
    mergeAndConvertToOsfls(pfss_oi[0], scs_oi[0]) # Sun-earth connection merge    
    shutil.copy(rootpath + wsafile[0], writepath + "/WSA_OUT") # this guy just copies the wsa-file
    count += 1
    print('Finished converting date {} after {} seconds.'.format(timestamp,time.time()-start_time))
    print('File {} out of {}'.format(count,nDates))      

print("Execution time: {} seconds".format(time.time()-start_time))