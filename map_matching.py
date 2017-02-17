import csv
import math
from math import cos, sin, asin, sqrt, atan2
from operator import attrgetter

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
<<<<<<< HEAD
        matchedLink = None
=======
        matchedLinke = None
    def __repr__(self):
        outstring = "\nSampleID: " + str(self.sampleID)
        outstring += "\nPoint3D: " + repr(self.point3D)
        outstring += "\nSpeed: " + str(self.speed)
        outstring += "\nHeading: " + str(self.heading)
        outstring += "\n"
        return outstring
>>>>>>> 102ec693545fd431f9ae83bd1e18a0c3f50e5c5f

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
        self.shapeInfo = parse_points_from_string(list_in[14], self)		# contains an array of shape entries consisting of the latitude and longitude (in decimal degrees) and elevation (in decimal meters) for the link's nodes and shape points ordered as reference node, shape points, non-reference node. The array entries are delimited by a vertical bar character and the latitude, longitude, and elevation values for each entry are delimited by a forward slash character (e.g. lat/lon/elev|lat/lon/elev). The elevation values will be null for links that don't have 3D data.
        self.curvatureInfo = list_in[15]		# contains an array of curvature entries consisting of the distance from reference node (in decimal meters) and curvature at that point (expressed as a decimal value of 1/radius in meters). The array entries are delimited by a vertical bar character and the distance from reference node and curvature values for each entry are separated by a forward slash character (dist/curvature|dist/curvature). This entire field will be null if there is no curvature data for the link.
        self.slopeInfo = list_in[16]		# contains an array of slope entries consisting of the distance from reference node (in decimal meters) and slope at that point (in decimal degrees). The array entries are delimited by a vertical bar character and the distance from reference node and slope values are separated by a forward slash character (dist/slope|dist/slope). This entire field will be null if there is no slope data for the link.
        self.refNode = self.shapeInfo[0]
        self.get_headings()

    def get_headings(self):
        index = 0
        for index in range(len(self.shapeInfo)):
            directions = []
            if index > 0:
                directions.append(azimuth(self.shapeInfo[index - 1], self.shapeInfo[index]))
            if index < len(self.shapeInfo) - 1:
                directions.append(azimuth(self.shapeInfo[index], self.shapeInfo[index+1]))
            final_direction = sum(directions)*1.0/len(directions)
            self.shapeInfo[index].heading = final_direction



def azimuth(point1, point2):
    y = sin(point2.longitude - point1.longitude)*cos(point2.latitude)
    x = cos(point1.latitude)*sin(point2.latitude) - sin(point1.latitude)*cos(point2.latitude)*cos(point2.longitude - point1.longitude)
    return 180*atan2(y, x)/math.pi

class Point3D:
    def __init__(self, latitude, longitude, elevation, parent=None):
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.parentRef = parent
        self.heading = None

    def distance_3D(self, other_point_3D):
        pass# not implemented yet

    def distance_2D(self, other):
        # caclulate 2D distance between two points
        # based on Haversine formula
        p = 0.017453292519943295 # pi / 180
        a = 0.5 - cos((other.latitude - self.latitude) * p)/2 + cos(self.latitude * p) * cos(other.latitude * p) * (1 - cos((other.longitude - self.longitude) * p)) / 2
        return 12742 * asin(sqrt(a))
    def __repr__(self):
        return  "(" + str(self.latitude) + ", " + str(self.longitude) + ")"
# parsing functions

def parse_points_from_string(shape_info, parent=None):
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
        return_obj = Point3D(latitude, longitude, elevation, parent)
        return_list.append(return_obj)
    return return_list

def parse_links_csv(path):
    print "Parsing Links from File"
    with open(path, 'rb') as links_file:
        line_reader = csv.reader(links_file)
        links_list = []
        for link_line in line_reader:
            links_list.append(Link(link_line))
    return links_list

def parse_nodes_csv(path, start, decimation, end):
    print "Parsing Nodes from File"
    probe_points_file = open(path, 'rb')
    probe_point_line_reader = csv.reader(probe_points_file)

    counter = start
    probe_list = []
    for probe_point_line in probe_point_line_reader:
        if counter % decimation != 0:
            counter = counter + 1
            continue
        # print counter
        counter = counter + 1
        # print probe_point_line
        probe_point = ProbePoint(probe_point_line)
        probe_list.append(probe_point)
        if counter >= end:
            break

    probe_points_file.close()
    return probe_list
# helper functions

def assign_nodes(nodepath, linkpath, probe_points, links):
    pass




def sort_control_points(points):
    # given a list of points of type Point3D, generates two lists sorted
    # by increasing latitude and longitude
    print "Sorting control points"
    sorted_latitudes = sorted(points, key=attrgetter('latitude'))
    sorted_longitudes = sorted(points, key=attrgetter('longitude'))
    return sorted_latitudes, sorted_longitudes

def write_links_csv(links, filename):
    # given a list of links of type Link, writes them to a CSV file of a given name
    out_file = open(filename, 'wb')
    writer = csv.writer(out_file, delimiter = ",")
    for link in links:
        for control_point in link.shapeInfo:
            line = [link.linkPVID, control_point.latitude, control_point.longitude]
            writer.writerow(line)
    out_file.close()

def write_points_csv(control_points, filename):
    # given a list of links of type Point3D, writes them to a CSV file of a given name
    # only writes the assigned linkPVID, latitude, and longitude
    out_file = open(filename, 'wb')
    writer = csv.writer(out_file, delimiter = ",")
    for control_point in control_points:
        line = [control_point.parentRef.linkPVID, control_point.latitude, control_point.longitude]
        writer.writerow(line)
    out_file.close()

def write_probe_csv(probe_points, filename):
    # given a list of links of type Point3D, writes them to a CSV file of a given name
    # only writes the assigned sampleID, latitude, and longitude
    out_file = open(filename, 'wb')
    writer = csv.writer(out_file, delimiter=",")
    for probe in probe_points:
        line = [probe.sampleID, probe.point3D.latitude, probe.point3D.longitude]
        writer.writerow(line)
    out_file.close()

def get_control_point_list(links):
    # given a list of type Link, returns all of the control points in those links
    points = []
    for link in links:
        points += link.shapeInfo
    return points

def control_points_in_range(probe_point, radius_meters, sorted_latitudes, sorted_longitudes):
    # given a specific probe point (Point3D), finds all of the control points (point3D) within a specific range
    # the range is given in meters, and defines a square (roughly speaking)
    # also requires the two lists of control points sorted by latitude and longitude (lists of Point3D)

    point3d = probe_point.point3D
    radius_latitude = meters_to_degrees_latitude(radius_meters)
    lat_start = point3d.latitude - radius_latitude
    lat_start_index = bisect_point3D(sorted_latitudes, lat_start, "latitude")
    lat_end = point3d.latitude + radius_latitude
    print "latitude range:", lat_start, lat_end
    lat_end_index = bisect_point3D(sorted_latitudes, lat_end, "latitude")
    print "latitude indices:", lat_start_index, lat_end_index
    radius_longitude = meters_to_degrees_longitude(radius_meters, point3d.latitude)
    long_start = point3d.longitude - radius_longitude
    long_start_index = bisect_point3D(sorted_longitudes, long_start, "longitude")
    long_end = point3d.longitude + radius_longitude
    print "longitude range:", long_start, long_end
    long_end_index = bisect_point3D(sorted_longitudes, long_end, "longitude")
    print "longitude indices:", long_start_index, long_end_index
    potential_long = sorted_latitudes[lat_start_index:lat_end_index]
    potential_lat = sorted_longitudes[long_start_index:long_end_index]
    potential_points = list(set(potential_lat) & set(potential_long))
    return potential_points

def bisect_point3D(points, value, attribute):
    # wrapper function for recursive call
    # given a list of point3D and a search value (either lat or long)
    # and the attribute indicating whether it's lat or long,
    # returns the index closest to the value via recursive binary search
    return bisect_point3D_helper(points, value, attribute, 0, len(points))

def bisect_point3D_helper(points, value, attribute, start_index, end_index):
    mid_index = (start_index + end_index)/2
    if end_index-start_index is 1 or end_index - start_index is 0:
        return start_index
    if attribute is "latitude":
        mid = points[mid_index].latitude
    elif attribute is "longitude":
        mid = points[mid_index].longitude
    if value > mid:
        return bisect_point3D_helper(points, value, attribute, mid_index+1, end_index)
    elif value < mid:
        return bisect_point3D_helper(points, value, attribute, start_index, mid_index)
    else:
        return mid_index

def heading_diff(probe_point, control_point):
    # aligns the probe point with a direction on the control point
    # computes the difference between the heading of the control point
    # and the probe point

    if control_point.parentRef.directionOfTravel is "F":
        headings = [control_point.heading]
        direction = "F"
    if control_point.parentRef.directionOfTravel is "T":
        headings = [180 + control_point.heading]
        direction = "T"
    else:
        headings = [control_point.heading, 180+control_point.heading]
    diffs = [compare_angles(probe_point.heading, heading) for heading in headings]
    if len(headings) is 2:
        if diffs[0] > diffs[1]:
            direction = "T"
        else:
            direction = "F"
    return direction, min(diffs)

def compare_angles(angle1, angle2):
    # finds the minimum difference between two angles
    diff = abs(angle1 - angle2) % 360.0
    if diff > 180:
        diff = 360-diff
    return diff

#untested for now
def snap_to_link(probe_point, link):
    best_index = 0
    best_distance = 100000000
    best_x = 0
    best_y = 0

    py = latitude_degrees_to_meters(probe_point.point3D.latitude)
    px = longitude_degrees_to_meters(probe_point.point3D.longitude, probe_point.point3D.latitude)

    for i in range(len(link.shapeInfo) - 1):
        point1 = link.shapeInfo[i]
        point2 = link.shapeInfo[i+1]
        y1 = latitude_degrees_to_meters(point1.latitude)
        y2 = latitude_degrees_to_meters(point2.latitude)
        x1 = longitude_degrees_to_meters(point1.longitude, point1.latitude)
        x2 = longitude_degrees_to_meters(point2.longitude, point2.latitude)
        if i == 0:
            refpx = x1
            refpy = y1

        vx = x2 - x1
        vy = y2 - y1
        toPx = px - x1
        toPy = py - y1

        project = (toPx * vx + toPy * vy) / (vx*vx + vy*vy)
        projectVx = vx*project
        projectVy = vy*project

        link_distance_sqr = vx*vx + vy*vy
        signed_distance_sqr = projectVx*projectVx + projectVy*projectVy
        #negative dot product means opposite direction
        if project < 0:
            signed_distance_sqr *= -1

        # add projected vector to the reference point.
        if signed_distance_sqr < 0:
            projectx = x1
            projecty = y1
        elif signed_distance_sqr > link_distance_sqr:
            projectx = x2
            projecty = y2
        else:
            projectx = x1 + projectVx
            projecty = y1 + projectVy

        #calculate whether this is better than the last
        distance = sqrt((px - projectx)**2 + (px - projecty)**2)
        if distance < best_distance:
            print "distsqr: ", signed_distance_sqr
            print "controldistsqr: ", link_distance_sqr
            best_index = i
            best_distance = distance
            best_x = projectx
            best_y = projecty
    probe_point.distFromLink = best_distance
    probe_point.distFromRef = sqrt((best_x - refpx)**2 + (best_y - refpy)**2)
    print "best point: ", best_x, ", ", best_y
    print "probe point: ", px, ", ", py
    print "control points: ", link.shapeInfo[best_index], ", ", link.shapeInfo[best_index+1]
    print "probe point lat/long: ", probe_point.point3D.latitude, ", ", probe_point.point3D.longitude
    best_lat = meters_to_degrees_latitude(best_y)
    best_long = meters_to_degrees_longitude(best_x, best_lat)
    return Point3D(best_lat, best_long, 0)

def speed_diff(probe_point, control_point, direction):
    # finds the difference between the speed of the probe_point (Point3D)
    # and the speed limit associated with the control_point (Point3D)
    # in the direction of travel for the probe_point
    # (requires the direction of travel returned by heading_score)

    if direction == "T":
        link_speed = control_point.parentRef.toRefSpeedLimit
    else:
        link_speed = control_point.parentRef.fromRefSpeedLimit
    return abs(probe_point.speed - link_speed)

def latitude_degrees_to_meters(latitude_distance):
    # converts latitude in degrees to meters
    return latitude_distance*111111

def meters_to_degrees_latitude(meters):
    # convert distance in meters to latitude in degrees
    return meters/111111.0

def longitude_degrees_to_meters(longitude_distance, current_latitude):
    # converts longitude in degrees to meters
    return longitude_distance*111111*cos(current_latitude)

def meters_to_degrees_longitude(meters, current_latitude):
    # convert distance in meters to longitude in degrees
    return meters/111111.0/cos(current_latitude)

def normalize(differences, max_diff):
    # given a list of differences, returns a list of scores
    # assuming that low difference is good and high difference is bad
    # score is between 0 and 1, 0 being good (low difference) and 1 being bad (high difference)
    return [1.0*diff/max_diff for diff in differences]

jacob_path = '/Users/jdbruce/Downloads/WQ2017/Geospatial/probe_data_map_matching/'
will_path = 'c:/Users/Will Molter/Documents/College/Winter 2017/EECS 395/proj2/'
links_filename = 'Partition6467LinkData.csv'
points_filename = 'Partition6467ProbePoints.csv'
path = jacob_path

links = parse_links_csv(path + links_filename)
control_points = get_control_point_list(links)
sorted_latitudes, sorted_longitudes = sort_control_points(control_points)

points = parse_nodes_csv(path + points_filename, 0, 1, 1)
# snapped_points = []
# for probe_point in points:
#     candidate_points = control_points_in_range(probe_point, 500, sorted_latitudes, sorted_longitudes)
#     best_link = candidate_points[0].parentRef
#     snapped = snap_to_link(probe_point, best_link)
#     snapped.parentRef = best_link
#     snapped_points.append(snapped)
#     probe_point.matchedLink = best_link

    # distances = [probe_point.point3D.distance_2D(control_point) for control_point in candidate_points]
    # speed_scores = [abs(probe_point.speed - control_point.parentRef.toRefSpeedLimit) for control_point in candidate_points]
    # heading_info = [heading_score(probe_point, control_point) for control_point in candidate_points]
    # heading_scores = [pair[1] for pair in heading_info]
    # directions = [pair[0] for pair in heading_info]
    # print "directions: ", directions
    # print "heading scores: ", heading_scores
    # best_link = None
    # probe_point.linkPVID = best_link.linkPVID
    # probe_point.matchedLink = best_link

points = parse_nodes_csv(path + points_filename, 0, 1, 2)

search_radius = 100

for probe_point in points:
    print probe_point
    candidate_points = control_points_in_range(probe_point, search_radius, sorted_latitudes, sorted_longitudes)
    print "Number of Candidate Points: ", len(candidate_points)

    # distance scoring
    distances = [probe_point.point3D.distance_2D(control_point) for control_point in candidate_points]
    print "Distances :", distances
    distance_scores = normalize(distances, search_radius)
    print "Distance Scores: ", distance_scores

    # heading scoring
    heading_info = [heading_diff(probe_point, control_point) for control_point in candidate_points]
    heading_diffs = [pair[1] for pair in heading_info]
    print "Heading Diffs: ", heading_diffs
    heading_scores = normalize(heading_diffs, 180)
    print "heading scores: ", heading_scores

    # speed scoring
    directions = [pair[0] for pair in heading_info]
    print "directions: ", directions
    speed_diffs = [speed_diff(probe_point, candidate_points[i], directions[i]) for i in range(len(directions))]
    print "speed_diffs", speed_diffs
    speed_scores = normalize(speed_diffs, 150)
    print "speed scores: ", speed_scores

    # score totaling
    distance_weight = 1.0
    heading_weight = 1.0
    speed_weight = 1.0

    scores = [distance_weight*distance_scores[i] + heading_weight*heading_scores[i] + speed_weight*speed_scores[i] for i in range(len(candidate_points))]
    print "Scores: ", scores

    best_score_index = scores.index(min(scores))
    print "Best index: ", best_score_index
    best_point = candidate_points[best_score_index]
    print "Best point: " + repr(best_point)

    # print "For sampleID: ", probe_point.sampleID
    # print "We matched link: ", probe_point.linkPVID
    #write_points_csv(candidate_points, "points.csv")



#probe_list = assign_nodes(will_path + points_filename, will_path + links_filename)

write_probe_csv(points, "probes.csv")
# write_points_csv(snapped_points, "snapped.csv")
write_links_csv([probe.matchedLink for probe in points], "links.csv")

# write_probe_csv(points, "probes.csv")
#write_links_csv([probe.matchedLink for probe in probe_list], "links.csv")


# parse_point_from_string('51.4965800/9.3862299/|51.4994700/9.3848799/')
