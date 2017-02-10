import csv
from math import cos, asin, sqrt

# classes

class ProbePoint:
    def __init__(self, list_in):
        self.sampleID = int(list_in[0]) 	# is a unique identifier for the set of probe points that were collected from a particular phone.
        self.dateTime = list_in[1]          # is the date and time that the probe point was collected.
        self.sourceCode = int(list_in[2])  	# is a unique identifier for the data supplier (13 = COMPANY).
        self.point3D = Point3D(float(list_in[3]), float(list_in[4]), float(list_in[5])) # lat, long, elevation
        self.speed = int(list_in[6])     # is the speed in KPH.
        self.heading = int(list_in[7])   # is the heading in degrees.
        linkPVID = None         # is the published versioned identifier for the link.
        direction = None        # is the direction the vehicle was travelling on thelink (F = from ref node, T = towards ref node).
        distFromRef = None      # is the distance from the reference node to the map-matched probe point location on the link in decimal meters.
        distFromLink = None       # is the perpendicular distance from the map-matched probe point location on the link to the probe point in decimal meters.
        matchedLinke = None

class Link:
    def __init__(self, list_in):
        self.linkPVID = int(list_in[0])		# is the published versioned identifier for the link.
        self.refNodeID = int(list_in[1])
        self.nrefNodeID = int(list_in[2])		# is the internal identifier for the links non-reference node.
        self.length = float(list_in[3])			# is the length of the link (in decimal meters).
        self.functionalClass = int(list_in[4])		# is the functional class for the link (1-5).
        self.directionOfTravel = list_in[5]	# is the allowed direction of travel for the link (F from ref node, T towards ref node, B - both)
        self.speedCategory = int(list_in[6])		# is the speed category for the link (1-8).
        self.fromRefSpeedLimit = int(list_in[7])	# is the speed limit for the link (in kph) in the direction of travel from the reference node.
        self.toRefSpeedLimit = int(list_in[8])		# is the speed limit for the link (in kph) in the direction of travel towards the reference node.
        self.fromRefNumLanes = int(list_in[9])		# is the number of lanes for the link in the direction of travel from the reference node.
        self.toRefNumLanes = int(list_in[10])		# is the number of lanes for the link in the direction of travel towards the reference node.
        self.multiDigitized = list_in[11]		# is a flag to indicate whether or not the link is multiply digitized (T is multiply digitized, F is singly digitized).
        self.urban = list_in[12]			# is a flag to indicate whether or not the link is in an urban area (T is in urban area, F is in rural area).
        self.timeZone = float(list_in[13])		# is the time zone offset (in decimal hours) from UTC.
        self.shapeInfo = parse_points_from_string(list_in[14])		# contains an array of shape entries consisting of the latitude and longitude (in decimal degrees) and elevation (in decimal meters) for the link's nodes and shape points ordered as reference node, shape points, non-reference node. The array entries are delimited by a vertical bar character and the latitude, longitude, and elevation values for each entry are delimited by a forward slash character (e.g. lat/lon/elev|lat/lon/elev). The elevation values will be null for links that don't have 3D data.
        self.curvatureInfo = list_in[15]		# contains an array of curvature entries consisting of the distance from reference node (in decimal meters) and curvature at that point (expressed as a decimal value of 1/radius in meters). The array entries are delimited by a vertical bar character and the distance from reference node and curvature values for each entry are separated by a forward slash character (dist/curvature|dist/curvature). This entire field will be null if there is no curvature data for the link.
        self.slopeInfo = list_in[16]		# contains an array of slope entries consisting of the distance from reference node (in decimal meters) and slope at that point (in decimal degrees). The array entries are delimited by a vertical bar character and the distance from reference node and slope values are separated by a forward slash character (dist/slope|dist/slope). This entire field will be null if there is no slope data for the link.
        self.refNode = self.shapeInfo[0]

class Point3D:
    def __init__(self, latitude, longitude, elevation):
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
    def distance_3D(self, other_point_3D):
        pass# not implemented yet
    def distance_2D(self, other):
        # caclulate 2D distance between two points
        # based on Haversine formula
        p = 0.017453292519943295 # pi / 180
        a = 0.5 - cos((other.latitude - self.latitude) * p)/2 + cos(self.latitude * p) * cos(other.latitude * p) * (1 - cos((other.longitude - self.longitude) * p)) / 2
        return 12742 * asin(sqrt(a))

# parsing functions

def parse_points_from_string(shape_info):
    string_list = shape_info.split('|')
    return_list = []
    for point in string_list:
        lle = point.split('/')
        latitude = float(lle[0])
        longitude = float(lle[1])
        if lle[2] == '':
            elevation = None
        else:
            elevation = float(lle[2])
        return_obj = Point3D(latitude, longitude, elevation)
        return_list.append(return_obj)
    return return_list

def parse_links_csv(path):
    with open(path, 'rb') as links_file:
        line_reader = csv.reader(links_file)
        links_list = []
        for link_line in line_reader:
            links_list.append(Link(link_line))
    return links_list

# helper functions

def assign_nodes(nodepath, linkpath):
    print "Parsing links"
    links_list = parse_links_csv(linkpath)
    print "Done parsing links"
    with open(nodepath, 'rb') as probe_points_file:
        probe_point_line_reader = csv.reader(probe_points_file)

        counter = 0
        decimation = 100
        counter_stop = 1500
        probe_list = []
        for probe_point_line in probe_point_line_reader:
            if counter % decimation != 0:
                counter = counter + 1
                continue
            # print probe_point_line
            print counter
            probe_point = ProbePoint(probe_point_line)
            probe_list.append(probe_point)
            best_distance = float("inf")
            for link in links_list:
                this_distance = probe_point.point3D.distance_2D(link.refNode)
                if this_distance < best_distance:
                    best_distance = this_distance
                    best_link = link
            probe_point.linkPVID = best_link.linkPVID
            probe_point.matchedLink = best_link
            print "For sampleID: ", probe_point.sampleID
            print "We matched link: ", probe_point.linkPVID

            counter = counter + 1
            if counter > counter_stop:
                break
    return probe_list

def write_links_csv(links, filename):
    out_file = open(filename, 'wb')
    writer = csv.writer(out_file, delimiter = ",")
    for link in links:
        for control_point in link.shapeInfo:
            line = [link.linkPVID, control_point.latitude, control_point.longitude]
            writer.writerow(line)
    out_file.close()

def write_probe_csv(probe_points, filename):
    out_file = open(filename, 'wb')
    writer = csv.writer(out_file, delimiter=",")
    for probe in probe_points:
        line = [probe.sampleID, probe.point3D.latitude, probe.point3D.longitude]
        writer.writerow(line)
    out_file.close()



jacob_path = '/Users/jdbruce/Downloads/WQ2017/Geospatial/probe_data_map_matching/'
will_path = 'c:/Users/Will Molter/Documents/College/Winter 2017/EECS 395/proj2/'
links_filename = 'Partition6467LinkData.csv'
points_filename = 'Partition6467ProbePoints.csv'

probe_list = assign_nodes(will_path + points_filename, will_path + links_filename)
write_probe_csv(probe_list, "probes.csv")
write_links_csv([probe.matchedLink for probe in probe_list], "links.csv")

# parse_point_from_string('51.4965800/9.3862299/|51.4994700/9.3848799/')
