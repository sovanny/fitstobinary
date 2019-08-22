import math, numpy, struct, ctypes, time, shutil, collections, random
from os import walk, makedirs
from enum import Enum
from astropy.io import fits
from astropy.time import Time

################## EXTRACT TIMESTAMPS AND SAVE FILEPATHS #########################

writepath = "./" # set output rootfolder
# set names of file types and output folder names
PFSS_IO = "PFSS_IO"
PFSS_OI = "PFSS_OI"
SCS_OI = "SCS_OI"
WSA_OUT = "WSA_OUT"
WSA_VEL = "WSA_VEL"

SUB_EARTH = "SUB_EARTH"

def populatefiles(folderpath, fileending, daterange, type):
    for (dirpath, dirnames, filenames) in walk(folderpath):
        for file in filenames:    
            if(file.endswith(fileending)):
                date = file[daterange[0]:daterange[1]]
                if(type == 'PFSS_IO'): fitsfilesdates[date] = []
                try: fitsfilesdates[date].append((folderpath + '/' + file, type))
                except: continue       
        break

# POPULATE THE FILEPATH VARIABLE
# This assumes there are multiple timesteps to be processed, 
# and that their divided into folders for each filetype.
# But the script works for a single timestep as well.
# Enter parameters: folderpath, fileending, indices for timestamp in filename, filetype
fitsfilesdates = {} 
rootpath = "/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/"
populatefiles(rootpath + "trace_pfss_intoout",'fits.gz' , (0,12), PFSS_IO) # PFSS_IO has to be populated first
populatefiles(rootpath + "trace_pfss_outtoin", 'fits.gz', (0,12), PFSS_OI)
populatefiles(rootpath + "trace_scs_outtoin", 'fits.gz' , (0,12), SCS_OI)
populatefiles(rootpath + "WSA_OUT"          , 'fits'    , (4,16), WSA_OUT)
populatefiles(rootpath + "WSA_VEL"          , 'fits'    , (4,16), WSA_VEL)

# DELETE TIMESTEPS THAT DON'T HAVE ALL THE 5 REQUIRED FILES
datestodelete = []
for date in fitsfilesdates:
    if(len(fitsfilesdates[date]) != 5):
        datestodelete.append(date)
for date in datestodelete:
    del fitsfilesdates[date]
print("Deleted {} dates because they were incomplete sets".format(len(datestodelete)))

nDates = len(fitsfilesdates)
print("Total {} dates to convert.".format(nDates))
ordereddates = collections.OrderedDict(sorted(fitsfilesdates.items())) # sort

# Example of what the Python Dict looks like for 3 timesteps
# For each entry (key: timestamp) the filepath and the type (tuples) for the 5 files are stored  
# {
#     '201708061158': [
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_pfss_intoout/201708061158R000_trace_pfss_intoout_tracing.fits.gz', 'PFSS_IO'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_pfss_outtoin/201708061158R000_trace_pfss_outtoin_tracing.fits.gz', 'PFSS_OI'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_scs_outtoin/201708061158R000_trace_scs_outtoin_tracing.fits.gz', 'SCS_OI'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/WSA_OUT/wsa_201708061158R000_gong.fits', 'WSA_OUT'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/WSA_VEL/vel_201708061158R000_gong.fits', 'WSA_VEL')
#         ], 
#     '201707150128': [
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_pfss_intoout/201707150128R000_trace_pfss_intoout_tracing.fits.gz', 'PFSS_IO'),
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_pfss_outtoin/201707150128R000_trace_pfss_outtoin_tracing.fits.gz', 'PFSS_OI'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_scs_outtoin/201707150128R000_trace_scs_outtoin_tracing.fits.gz', 'SCS_OI'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/WSA_OUT/wsa_201707150128R000_gong.fits', 'WSA_OUT'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/WSA_VEL/vel_201707150128R000_gong.fits', 'WSA_VEL')
#         ], 
#     '201707300236': [
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_pfss_intoout/201707300236R000_trace_pfss_intoout_tracing.fits.gz', 'PFSS_IO'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_pfss_outtoin/201707300236R000_trace_pfss_outtoin_tracing.fits.gz', 'PFSS_OI'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/trace_scs_outtoin/201707300236R000_trace_scs_outtoin_tracing.fits.gz', 'SCS_OI'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/WSA_OUT/wsa_201707300236R000_gong.fits', 'WSA_OUT'), 
#         ('/Volumes/OPENSPACE/WSA_OUT_FIELDLINES_SEPMOD/WSA_VEL/vel_201707300236R000_gong.fits', 'WSA_VEL')
#         ]
# }

################### FUNCTION DEFINITIONS for processing fieldlines ###################

def sph2cart(coord):
    return [
         coord[2] * math.sin(coord[1]) * math.cos(coord[0]), 
         coord[2] * math.sin(coord[1]) * math.sin(coord[0]),
         coord[2] * math.cos(coord[1])
    ]

def coord2datarow(coord):
    sph = sph2cart(coord[0:3]) # spherical coords
    sph.append(coord[5]) # magnitude for field strength, wiht polarity
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

# Find closed fieldlines, based on differing polarities in first and last point of fieldline
# Returns the indices of the fieldlines that are closed
def pickClosed(filepath):
    fl_fits = fits.open(filepath)
    fl_data = fl_fits[0].data
    fl_fits.close()
    indices_to_save = []

    for i in range(16200):
        b_radii = [pt[5] for pt in fl_data[i] if pt[0] > -900] 
        if (len(b_radii) < 2): continue
        first_b_radius = b_radii[0]
        last_b_radius = b_radii[-1]
        if(first_b_radius*last_b_radius < 0): #if product is negative/opposite signs
            indices_to_save.append(i)
    return indices_to_save

# Pick out field lines based on strength of field in first point of every line 
def pick_strength_based(filepath):
    fl_fits = fits.open(filepath)
    fl_data = fl_fits[0].data
    fl_fits.close()
    abs_bfieldvalues = [numpy.abs(line[0][5]) for line in fl_data[:16200]]
    del fl_data

    zscore = 2              # n of standard deviations
    hard_threshold = 15     # gauss
    n_extra_lines = 300     # lines to pick besides from the ones above threshold 
    # set the threshold to z standard deviations away from mean
    std_threshold = numpy.mean(abs_bfieldvalues) + zscore * numpy.std(abs_bfieldvalues)
    final_threshold = max(hard_threshold, std_threshold) # choose the hard threshold if std is low

    distribution_array = []
    indices_to_keep = []
    for i in range(16200):
        strength = abs_bfieldvalues[i]
        if(strength > final_threshold): indices_to_keep.append(i)
        # If not above limit, replicate element a number of times proportionate to field strength
        # Will be used to do random picking based on distribution later (Monte Carlo)
        else: distribution_array.extend( [i]*math.ceil(strength) ) 
    # print("Number of strong values for " + filepath + " is " + str(len(indices_to_keep)))

    # Pick an additonal set of lines based on distribution. If there is low activity, the above 
    # algorithm will probably result in few or zero lines, so this fills out with lines.
    count = 0
    while(count != n_extra_lines):
        random_index = random.choice(distribution_array)
        if random_index in indices_to_keep: continue
        indices_to_keep.append(random_index)
        count += 1

    return [indices_to_keep, 'mntcrl250']

# Pick lines where two polarities meet and make up the current sheet
def pick_sheet_lines(wsa_file):
    wsa_fits = fits.open(wsa_file)
    wsa_data = wsa_fits[0].data
    wsa_data = wsa_data[0] # the first layer has the coronal field information
    wsa_fits.close()

    sheet_indices = []
    for i in range(90):
        for j in range(180):
            if(i != 89 and numpy.sign(wsa_data[i][j]) != numpy.sign(wsa_data[i+1][j])):
                index = i*180  +j
                sheet_indices.append(index)
                index = (i+1)*180 +j
                sheet_indices.append(index)
            if(j != 179) and numpy.sign(wsa_data[i][j]) != numpy.sign(wsa_data[i][j+1]):
                index = i*180  +j
                sheet_indices.append(index)
                index = i*180  +j+1
                sheet_indices.append(index)

    return list(numpy.unique(sheet_indices))

# Pick the lines where closed and open field lines meet, so that the lines picked create a boundary.
# Do this by comparing every line end point in the pfss out-to-in set to every start and end point of 
# the closed lines in the pfss in-to-out set. If they are close to eachother coordinate wise, they 
# are classed as a boundary line.
def boundaryLines(filepath_io, filepath_oi, indices_closed):
    fl_io = (fits.open(filepath_io))[0].data
    fl_oi = (fits.open(filepath_oi))[0].data

    threshold = math.sqrt(2)
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
        # get the last coordinates for out-to-in (the ones nearest the sun)
        oi_last_phi = math.degrees(fl_oi[i][last][0])
        oi_last_theta = math.degrees(fl_oi[i][last][1])
        for j in indices_closed:
            # calculate distances and compare
            dist_last_phi = abs(io_last_phi[j] - oi_last_phi)
            dist_last_theta = abs(io_last_theta[j] - oi_last_theta)
            dist_first_phi = abs(io_first_phi[j] - oi_last_phi)
            dist_first_theta = abs(io_first_theta[j] - oi_last_theta)
            if(dist_last_phi < threshold and dist_last_theta < threshold):
                boundary_lines_oi.append(i)
                break
            if(dist_first_phi < threshold and dist_first_theta < threshold):
                boundary_lines_oi.append(i)
                break
    return [boundary_lines_oi, 'boundary']
# end def pickBoundaryLines()

# Pick every nth between the selected boundary indices to fill out the field line set
def fillOut(boundary_indices, step):
    output = []
    previous_index = 0
    for index in boundary_indices:
        subrange = list(range(previous_index, index, step))
        output.extend(subrange)
        previous_index = index
    subrange = list(range(index, 16200, step))
    output.extend(subrange)
    return [output, 'boundary_filled']

# Produce the .osfls (OpenSpace Field Line Sequence) binary
def toOsfls(filepath, modelname, indices, typename, extra_indices, extra_file):
    fl_fits = fits.open(filepath)
    fl_data = fl_fits[0].data
    fl_fits.close()
    
    if(modelname != 'PFSS_IO'):
        vel_fits = fits.open(extra_file)
        vel_data = vel_fits[0].data
        vel_data = vel_data[1] # the second layer has the wind speed
        vel_fits.close()
        indices.extend(extra_indices)
        indices = list(numpy.unique(indices))
        
    versionNumber = 0
    [triggerTime, pathSafeTimeString] = obstimeToJ2000(fl_fits[0].header['OBSTIME'])
    filepath = pathSafeTimeString + '.osfls'
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
            if(i in extra_indices): openness = 0
            elif(points[0][3] < 0): openness = 1 # open and negative
            else: openness = 2 # open and positive
            xtra = [openness]*len(points) # same value for all points
            extraQuantities_1.extend(xtra)
        else:
            velocity = vel_data[math.floor(i/180)][i%180] # speed at 21.5
            if(points[0][3] < 0): velocity = (-1)*velocity # negative polarity
            xtra = [velocity]*len(points) # same value for all points
            extraQuantities_1.extend(xtra)
            in_sheet = 0 
            if(i in extra_indices): in_sheet = 1 # 1 if in current sheet 
            xtra = [in_sheet]*len(points) # same value for all points
            extraQuantities_2.extend(xtra)
    
    nLines = len(lineStart)
    nExtras = 1
    extraQuantities = extraQuantities_1
   
    if(modelname == 'PFSS_IO'):
        extraQuantityNames = ['Open/closed lines \0']
    else: 
        nExtras = 2
        extraQuantities.extend(extraQuantities_2)
        extraQuantityNames = ['Polarity and solar wind speed \0', 'Current sheet \0']
    
    nStringBytes = sum([len(s) for s in extraQuantityNames])
    allNamesInOne = ''
    for s in extraQuantityNames:
        allNamesInOne += s
    
    # Prepare data for writing to binary file. Using Struct and pack
    typestr = '= i d i ? Q Q Q Q %sl %sL %sf %sf %ss' % (nLines, nLines, 3*nVert, nExtras*nVert, nStringBytes)
    struct_to_write = struct.Struct(typestr)
    values_to_write = (versionNumber, triggerTime, model, isMorphable, nLines, nVert, nExtras, nStringBytes)
    values_to_write += (*lineStart, *lineCount, *vertexPositions, *extraQuantities, allNamesInOne.encode('utf-8'))
    buffer = ctypes.create_string_buffer(struct_to_write.size)  
    struct_to_write.pack_into(buffer, 0, *values_to_write)
    fout = open(writepath + modelname + '/' + filepath, 'wb')
    fout.write(buffer)
    fout.close()
# end def toOsfls()

# Produce the .osfls (OpenSpace Field Line Sequence) binary
# - specially for the sub-satellite points, aka sub-earth points, data set . 
# Here, the SCS and PFSS out-to-in sets are merged for the last 540 indices. 
def mergeAndConvertToOsfls(filepathPFSS, filepathSCS):
    fl_fits_pfss = fits.open(filepathPFSS)
    fl_data_pfss = fl_fits_pfss[0].data
    fl_fits_pfss.close()
    fl_fits_scs = fits.open(filepathSCS)
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
    filepath = pathSafeTimeString + '.osfls'
    model = Model.Wsa.value
    isMorphable = False
    
    nVert = 0
    lineStart = []
    lineCount = []
    vertexPositions = []
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
  
        polarity = 0
        if(points_scs[0][3] >= 0):
            polarity = 1
        xtra = [polarity]*combinedVertLen
        extraQuantities_polarity.extend(xtra)
       
        subsatellitepoint = 1
        if (i in range(16380, 16560)):
            subsatellitepoint = 0
        elif (i in range(16560, 16740)):
            subsatellitepoint = 2
        xtra = [subsatellitepoint]*combinedVertLen
        extraQuantities_level.extend(xtra)
       
    nLines = len(lineStart)

    nExtras = 2
    extraQuantities = extraQuantities_polarity
    extraQuantities.extend(extraQuantities_level)

    extraQuantityNames = ['Polarity \0', 'Sub-Earth points \0']

    nStringBytes = sum([len(s) for s in extraQuantityNames])
    allNamesInOne = ''
    for s in extraQuantityNames:
        allNamesInOne += s
    
    # Prepare data for writing to binary. Using Struct and pack
    typestr = '= i d i ? Q Q Q Q %sl %sL %sf %sf %ss' % (nLines, nLines, 3*nVert, nExtras*nVert, nStringBytes)
    struct_to_write = struct.Struct(typestr)
    values_to_write = (versionNumber, triggerTime, model, isMorphable, nLines, nVert, nExtras, nStringBytes)
    values_to_write += (*lineStart, *lineCount, *vertexPositions, *extraQuantities, allNamesInOne.encode('utf-8'))
    
    buffer = ctypes.create_string_buffer(struct_to_write.size)    
    struct_to_write.pack_into(buffer, 0, *values_to_write)
    
    fout = open(writepath + SUB_EARTH + '/' + filepath, 'wb')
    fout.write(buffer)
    fout.close()
# end def mergeAndConvertToOsfls()

###################### PROCESS ALL FILES ##############################
start_time = time.time()
count = 0

# Create the output directories
makedirs(writepath + PFSS_IO, exist_ok=True)
makedirs(writepath + PFSS_OI, exist_ok=True)
makedirs(writepath + SCS_OI, exist_ok=True)
makedirs(writepath + WSA_OUT, exist_ok=True)
makedirs(writepath + SUB_EARTH, exist_ok=True)

# relies heavily on having a sorted list of files for each timestep
# so that pfss_io goes first and picks out the lines
for timestamp in ordereddates:
    pfss_io = fitsfilesdates[timestamp][0]
    pfss_oi = fitsfilesdates[timestamp][1]
    scs_oi = fitsfilesdates[timestamp][2]
    wsafile = fitsfilesdates[timestamp][3]
    velfile = fitsfilesdates[timestamp][4]
    indices_closed_lines = pickClosed(pfss_io[0])
    indices = pick_strength_based(pfss_io[0])
    toOsfls(pfss_io[0], pfss_io[1], indices[0], indices[1], indices_closed_lines, 0)
    print('Finished converting PFSS in-to-out after {} seconds.'.format(time.time()-start_time))
    indices_boundary = boundaryLines(pfss_io[0], pfss_oi[0], indices_closed_lines)
    print('Picked boundary lines after after {} seconds.'.format(time.time()-start_time))
    indices = fillOut(indices_boundary[0], 23) # change step here for sparser
    sheet_indices = pick_sheet_lines(wsafile[0])
    print('Picked current sheet lines after after {} seconds.'.format(time.time()-start_time))
    toOsfls(pfss_oi[0], pfss_oi[1], indices[0], indices[1], sheet_indices, velfile[0])
    print('Finished converting PFSS out-to-in after {} seconds.'.format(time.time()-start_time))
    toOsfls(scs_oi[0], scs_oi[1], indices[0], indices[1], sheet_indices, velfile[0])
    print('Finished converting SCS out-to-in after {} seconds.'.format(time.time()-start_time))
    mergeAndConvertToOsfls(pfss_oi[0], scs_oi[0]) # Sub-earth track merge    
    shutil.copy(wsafile[0], writepath + "/WSA_OUT") # this line just copies the wsa-file
    count += 1
    print('Finished converting date {} after {} seconds.'.format(timestamp,time.time()-start_time))
    print('File {} out of {}'.format(count,nDates))      

print("Execution time: {} seconds".format(time.time()-start_time))