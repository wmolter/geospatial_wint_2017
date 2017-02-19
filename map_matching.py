import csv
import math
import numpy as np
from math import cos, sin, asin, sqrt, atan2
from operator import attrgetter

# classes

class ProbePoint:
    def __init__(self, list_in, ID=-1):
        self.uniqueID = ID
        self.sampleID = int(list_in[0]) 	# is a unique identifier for the set of probe points that were collected from a particular phone.
        self.dateTime = list_in[1]          # is the date and time that the probe point was collected.
        self.sourceCode = int(list_in[2])  	# is a unique identifier for the data supplier (13 = COMPANY).
        self.point3D = Point3D(float(list_in[3]), float(list_in[4]), float(list_in[5])) # lat, long, elevation
        self.speed = int(list_in[6])     # is the speed in KPH.
        self.heading = int(list_in[7])   # is the heading in degrees.
        if len(list_in) > 8 and list_in[8] is not "":
            self.linkPVID = int(list_in[8])         # is the published versioned identifier for the link.
            self.direction = list_in[9]        # is the direction the vehicle was travelling on thelink (F = from ref node, T = towards ref node).
            # THIS WILL BE AN INDEX BEFORE YOU CHANGE IT
            self.snappedControlPoint = int(list_in[10])
            self.distFromRef = float(list_in[11])      # is the distance from the reference node to the map-matched probe point location on the link in decimal meters.
            self.distFromLink = float(list_in[12])       # is the perpendicular distance from the map-matched probe point location on the link to the probe point in decimal meters.
            self.uniqueID = int(list_in[13])
        else:
            self.linkPVID = None         # is the published versioned identifier for the link.
            self.direction = None        # is the direction the vehicle was travelling on thelink (F = from ref node, T = towards ref node).
            self.snappedControlPoint = None
            self.distFromRef = None      # is the distance from the reference node to the map-matched probe point location on the link in decimal meters.
            self.distFromLink = None       # is the perpendicular distance from the map-matched probe point location on the link to the probe point in decimal meters.
        self.matchedLink = None

    def __repr__(self):
        outstring = "\nSampleID: " + str(self.sampleID)
        outstring += "\nPoint3D: " + repr(self.point3D)
        outstring += "\nSpeed: " + str(self.speed)
        outstring += "\nHeading: " + str(self.heading)
        outstring += "\n"
        return outstring

    def to_csv_list(self):
        # if self.matchedLink is not None:
        #     # control_index = self.matchedLink.shapeInfo.index(self.snappedControlPoint)
        # else:
        control_index = -1
        return [self.sampleID, self.dateTime, self.sourceCode, self.point3D.latitude, self.point3D.longitude, self.point3D.elevation,
                self.speed, self.heading, self.linkPVID, self.direction, control_index,self.distFromRef, self.distFromLink, self.uniqueID]

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
        if list_in[16] is None or list_in[16] is "":
            self.slopeInfo = None
        else:
            self.slopeInfo = list_in[16].split("|")	# contains an array of slope entries consisting of the distance from reference node (in decimal meters) and slope at that point (in decimal degrees). The array entries are delimited by a vertical bar character and the distance from reference node and slope values are separated by a forward slash character (dist/slope|dist/slope). This entire field will be null if there is no slope data for the link.
        # contains an array of slope entries consisting of the distance from reference node (in decimal meters) and slope at that point (in decimal degrees). The array entries are delimited by a vertical bar character and the distance from reference node and slope values are separated by a forward slash character (dist/slope|dist/slope). This entire field will be null if there is no slope data for the link.
        self.elevSum = [0] * len(self.shapeInfo) # contains an array of sums of elevation, each corresponding to the control point
        self.elevCount = [0] * len(self.shapeInfo) # contains an array of counts of elevation data, each corresponding to a control point
        self.estSlopeInfo = None
        self.refNode = self.shapeInfo[0]
        # self.matchedProbePoints = []
        self.get_headings()

    def get_error(self):
        if self.slopeInfo is None or self.estSlopeInfo is None:
            return None
        total = 0
        for i in range(len(self.estSlopeInfo)):
            dist, slope = self.slopeInfo[i].split("/")
            slope = float(slope)
            estSlope = self.estSlopeInfo[i]
            if estSlope is None:
                print "slope unassigned"
            total += abs(estSlope-slope)
        average = total / len(self.estSlopeInfo)
        return average



    # def to_csv_list(self):
    #     return [self.linkPVID, [point.uniqueID for point in self.matchedProbePoints]]

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
        return

    def estimate_slope_info(self):
        # assumes that at least two points have already been snapped to this link, i.e. len(matchedProbePoints > 1)
        # assumes we can't use elevation of control points to calculate slope, only snapped probe points
        # assumes that elevSum and elevCount are both populated, i.e. all the control points have been run

        # first, get a list of control point dicts
        # each dict has points (a list), elevation, slope, and control point data

        # check to make sure that there are at least two non-zero entries in the elevation
        # data. Without at least that much, no slope information can be used

        count_count = 0
        for count in self.elevCount:
            if count > 0:
                count_count += 1
        if count_count < 2:
            return # self.estSlopeInfo is already none

        control_point_data = []

        for (control_point_index, control_point) in enumerate(self.shapeInfo):
            datum = {}
            if self.elevCount[control_point_index] is 0:
                datum["elevation"] = None
            else:
                datum["elevation"] = 1.0 * self.elevSum[control_point_index] / self.elevCount[control_point_index]
            datum["slope"] = None
            datum["control point"] = control_point
            control_point_data.append(datum)

        # estimate slopes between points by drawing lines between points with estimated elevations

        index = 0
        slopes = []

        current_slope = {"start_index": None, "end_index": None} # structure is: start index, end index, slope
        start_elevation = None

        for (point_index, datum) in enumerate(control_point_data):
            if datum["elevation"] is None:
                continue
            # elevation is defined for this point,
            # we should count it for calculating slope
            if current_slope["start_index"] is None:
                # this is the first point
                start_elevation = datum["elevation"]
                current_slope["start_index"] = point_index
            else:
                # this is the second point of one line, and the first of the next
                current_slope["end_index"] = point_index
                current_slope["slope"] = calculate_slope(self.shapeInfo[current_slope["start_index"]], start_elevation, self.shapeInfo[point_index], datum["elevation"])
                if current_slope["slope"] is None:
                    print "calculated slope is None"
                slopes.append(current_slope)
                # start the next line
                current_slope = {}
                current_slope["start_index"] = point_index


        # estimate slopes at points by taking the average of forward and backward slope, if available
        # if not available, take the average of any slope passing through the point

        for point_index in range(len(control_point_data)):
            # find whether it occurs at all in the slope list
            probe_point_datum = control_point_data[point_index]
            for (slope_index, slope_datum) in enumerate(slopes):
                if point_index <= slope_datum["start_index"]:
                    # the point is less than the first slope
                    # or it would have been handled already
                    # so estimate its slope by the first slope that we can find
                    probe_point_datum["slope"] = slope_datum["slope"]
                    break
                if (point_index > slope_datum["start_index"]) and (point_index < slope_datum["end_index"]):
                    # if it's between two indices, it's on the line (it didn't have elevation data)
                    # so assign it the slope of that line
                    probe_point_datum["slope"] = slope_datum["slope"]
                    break
                if point_index is slope_datum["end_index"]:
                    # this point is the end of one line
                    try:
                        assert(point_index is slopes[slope_index + 1]["start_index"])
                        # if there is another slope entry, this point should be the start of that as well
                        # this should be the most common case, where we average between the forward and backward slopes
                        probe_point_datum["slope"] = 1.0* (slope_datum["slope"] + slopes[slope_index + 1]["slope"]) / 2
                        break
                    except IndexError:
                        # if there was an Index error, this slope was actually the last one
                        # in that case, the slope should just be the backward slope
                        probe_point_datum["slope"] = slope_datum["slope"]
                        break
                if point_index > slope_datum["end_index"]:
                    # the point is further than the last slope entry we have
                    # so we estimate its slope by the last slope
                    probe_point_datum["slope"] = slope_datum["slope"]
                    break
        for datum in control_point_data:
            if datum["slope"] is None:
                print datum
                print "LinkID", self.linkPVID
        self.estSlopeInfo = [datum["slope"] for datum in control_point_data]

        return


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
        return 1000 * 12742 * asin(sqrt(a))
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

    counter = 0
    probe_list = []
    for probe_point_line in probe_point_line_reader:
        if counter % decimation != 0 or counter < start:
            counter = counter + 1
            continue
        # print counter
        counter = counter + 1
        # print probe_point_line
        probe_point = ProbePoint(probe_point_line, counter)
        probe_list.append(probe_point)
        if counter >= end and end is not -1:
            break

    probe_points_file.close()
    return probe_list
# helper functions

def calculate_slope(pt3D1, elev1, pt3D2, elev2):
    # calculates the slope as pt2 - pt1
    # if the slope goes up from pt1 to pt2, the value will be positive


    delta_elev = elev2 - elev1
    crow_distance = pt3D1.distance_2D(pt3D2) # 2D distance between the points

    output = math.degrees(math.atan2(delta_elev, crow_distance))
    if output > 10:
        print "Large slope: ", output
        print "estimated elevations: ", elev1, elev2
        print "actual elevations: ", pt3D1.elevation, pt3D2.elevation
        print "Points: ", pt3D1, pt3D2
        print "2D Distance: ", crow_distance
    return output

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
        line = [control_point.parentRef.linkPVID, control_point.latitude, control_point.longitude, control_point.heading]
        writer.writerow(line)
    out_file.close()

def write_probe_csv(probe_points, filename):
    # given a list of links of type Point3D, writes them to a CSV file of a given name
    # only writes the assigned sampleID, latitude, and longitude
    out_file = open(filename, 'wb')
    writer = csv.writer(out_file, delimiter=",")
    for probe in probe_points:
        line = [probe.sampleID, probe.point3D.latitude, probe.point3D.longitude, probe.heading]
        writer.writerow(line)
    out_file.close()

def write_probe_info_csv(probe_points, filename):
    out_file = open(filename, 'wb')
    writer = csv.writer(out_file, delimiter=",")
    for probe in probe_points:
        writer.writerow(probe.to_csv_list())
    out_file.close()

# def write_link_info_csv(links, filename):
#     out_file = open(filename, 'wb')
#     writer = csv.writer(out_file, delimiter=",")
#     for link in links:
#         writer.writerow(link.to_csv_list())
#     out_file.close()

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
    # print "latitude range:", lat_start, lat_end
    lat_end_index = bisect_point3D(sorted_latitudes, lat_end, "latitude")
    # print "latitude indices:", lat_start_index, lat_end_index
    radius_longitude = meters_to_degrees_longitude(radius_meters, point3d.latitude)
    long_start = point3d.longitude - radius_longitude
    long_start_index = bisect_point3D(sorted_longitudes, long_start, "longitude")
    long_end = point3d.longitude + radius_longitude
    # print "longitude range:", long_start, long_end
    long_end_index = bisect_point3D(sorted_longitudes, long_end, "longitude")
    # print "longitude indices:", long_start_index, long_end_index
    potential_lat = sorted_latitudes[lat_start_index:lat_end_index]
    potential_long = sorted_longitudes[long_start_index:long_end_index]
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

def azimuth(point1, point2):
    y = sin(point2.longitude - point1.longitude)*cos(point2.latitude)
    x = cos(point1.latitude)*sin(point2.latitude) - sin(point1.latitude)*cos(point2.latitude)*cos(point2.longitude - point1.longitude)
    return 180*atan2(y, x)/math.pi

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

    py = probe_point.point3D.latitude
    px = probe_point.point3D.longitude

    for i in range(len(link.shapeInfo) - 1):
        point1 = link.shapeInfo[i]
        point2 = link.shapeInfo[i+1]
        y1 = point1.latitude
        y2 = point2.latitude
        x1 = point1.longitude
        x2 = point2.longitude
        if i == 0:
            refpx = x1
            refpy = y1

        vx = longitude_degrees_to_meters(x2 - x1, y1)
        vy = latitude_degrees_to_meters(y2 - y1)
        toPx = longitude_degrees_to_meters(px - x1, y1)
        toPy = latitude_degrees_to_meters(py - y1)

        link_distance_sqr = vx*vx + vy*vy

        project = (toPx * vx + toPy * vy) / link_distance_sqr
        projectVx = vx*project
        projectVy = vy*project


        signed_distance_sqr = projectVx*projectVx + projectVy*projectVy
        #negative dot product means opposite direction

        # add projected vector to the reference point.
        if project < 0:
            signed_distance_sqr *= -1
            projectx = x1
            projecty = y1
        elif signed_distance_sqr > link_distance_sqr:
            projectx = x2
            projecty = y2
        else:
            projectx = x1 + meters_to_degrees_longitude(projectVx, y1)
            projecty = y1 + meters_to_degrees_latitude(projectVy)

        #calculate whether this is better than the last
        xdist = longitude_degrees_to_meters(px - projectx, y1)
        ydist = latitude_degrees_to_meters(py - projecty)
        distance = sqrt(xdist**2 + ydist**2)
        if distance < best_distance:
            # print "distsqr: ", signed_distance_sqr
            # print "controldistsqr: ", link_distance_sqr
            best_index = i
            best_distance = distance
            best_x = projectx
            best_y = projecty
    probe_point.distFromLink = best_distance
    probe_point.distFromRef = sqrt(longitude_degrees_to_meters(best_x - refpx, refpy)**2 + latitude_degrees_to_meters(best_y - refpy)**2)
    # print "best point: ", best_x, ", ", best_y
    # print "probe point: ", px, ", ", py
    # print "control points: ", link.shapeInfo[best_index], ", ", link.shapeInfo[best_index+1]
    # print "probe point lat/long: ", probe_point.point3D.latitude, ", ", probe_point.point3D.longitude
    return Point3D(best_y, best_x, 0)

def speed_diff(probe_point, control_point, direction):
    # finds the difference between the speed of the probe_point (Point3D)
    # and the speed limit associated with the control_point (Point3D)
    # in the direction of travel for the probe_point
    # (requires the direction of travel returned by heading_score)

    if direction == "T":
        link_speed = control_point.parentRef.toRefSpeedLimit
    else:
        link_speed = control_point.parentRef.fromRefSpeedLimit
    if link_speed is 999 or link_speed is 998:
        return 0
    return abs(probe_point.speed - link_speed)

def latitude_degrees_to_meters(latitude_distance):
    # converts latitude in degrees to meters
    return latitude_distance*111111

def meters_to_degrees_latitude(meters):
    # convert distance in meters to latitude in degrees
    return meters/111111.0

def longitude_degrees_to_meters(longitude_distance, current_latitude):
    # converts longitude in degrees to meters
    return longitude_distance*111111*cos(current_latitude*math.pi/180)

def meters_to_degrees_longitude(meters, current_latitude):
    # convert distance in meters to longitude in degrees
    return meters/111111.0/cos(current_latitude*math.pi/180)

def normalize(differences, max_diff):
    # given a list of differences, returns a list of scores
    # assuming that low difference is good and high difference is bad
    # score is between 0 and 1, 0 being good (low difference) and 1 being bad (high difference)
    return [1.0*diff/max_diff for diff in differences]

def match_probe_linkID(probe, link_dict):
    if probe.linkPVID is not None:
        probe.matchedLink = link_dict[probe.linkPVID]
        # probe.snappedControlPoint = probe.matchedLink.shapeInfo[probe.snappedControlPoint]
        # probe.matchedLink.matchedProbePoints.append(probe)
        probe.matchedLink.elevSum[probe.snappedControlPoint] += probe.point3D.elevation
        probe.matchedLink.elevCount[probe.snappedControlPoint] += 1

def get_links_dict(links_list):
    table = {}
    for link in links_list:
        table[link.linkPVID] = link
    return table

def get_total_error(links):
    total = 0
    max_error = 0
    max_link = 0
    errors = []
    for link in links:
        error = link.get_error()
        if error is not None:
            errors.append(error)
            if error > max_error:
                max_error = error
                max_link = link.linkPVID
    mean = sum(errors) / len(errors)
    median = np.median(np.array(errors))
    print "Number of links: ", len(errors)
    return mean, median, max_error, max_link

jacob_path = '/Users/jdbruce/Downloads/WQ2017/Geospatial/probe_data_map_matching/'
will_path = 'c:/Users/Will Molter/Documents/College/Winter 2017/EECS 395/proj2/'
links_filename = 'Partition6467LinkData.csv'
points_filename = 'Partition6467ProbePoints.csv'
points_info_file = "probe_info.csv"
path = will_path


links = parse_links_csv(path + links_filename)
control_points = get_control_point_list(links)
sorted_latitudes, sorted_longitudes = sort_control_points(control_points)


# points = parse_nodes_csv(path + points_filename, 100, 1, 101)
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

# points = parse_nodes_csv(path + points_filename, 250, 1, 251)
#
# search_radius = 200
#

#points = parse_nodes_csv(path + points_filename, 0, 1, -1)

search_radius = 500

print "Matching Points"
counter = 0
#out_file = open("probe_info.csv", 'wb')
#writer = csv.writer(out_file, delimiter=",")
in_file = open(path+points_filename, 'rb')
# in_file = open("probe_info.csv", 'rb')
probe_point_line_reader = csv.reader(in_file)

# for probe_point in points:
for probe_point_line in probe_point_line_reader:
    if counter % 1000 is 0:
        print counter
    # # print probe_point_line
    probe_point = ProbePoint(probe_point_line, counter)
    # print probe_point
    candidate_points = control_points_in_range(probe_point, search_radius, sorted_latitudes, sorted_longitudes)
    # all_candidates += candidate_points
    # print "Number of Candidate Points: ", len(candidate_points)

    # distance scoring
    distances = [probe_point.point3D.distance_2D(control_point) for control_point in candidate_points]
    # print "Distances :", distances
    distance_scores = normalize(distances, search_radius)
    # print "Distance Scores: ", distance_scores

    # heading scoring
    heading_info = [heading_diff(probe_point, control_point) for control_point in candidate_points]
    heading_diffs = [pair[1] for pair in heading_info]
    # print "Heading Diffs: ", heading_diffs
    heading_scores = normalize(heading_diffs, 180)
    # print "heading scores: ", heading_scores

    # speed scoring
    directions = [pair[0] for pair in heading_info]
    # print "directions: ", directions
    speed_diffs = [speed_diff(probe_point, candidate_points[i], directions[i]) for i in range(len(directions))]
    # print "speed_diffs", speed_diffs
    speed_scores = normalize(speed_diffs, 150)
    # print "speed scores: ", speed_scores

    # score totaling
    distance_weight = 5.0
    heading_weight = 1.0
    speed_weight = 1.0

    scores = [distance_weight*distance_scores[i] + heading_weight*heading_scores[i] + speed_weight*speed_scores[i] for i in range(len(candidate_points))]
    # print "Scores: ", scores

    if len(scores) is not 0:
        best_score_index = scores.index(min(scores))
        # print "Best index: ", best_score_index
        best_point = candidate_points[best_score_index]
        probe_point.matchedLink = best_point.parentRef
        # print "Best point: " + repr(best_point)

        matched_link = best_point.parentRef
        snapped = snap_to_link(probe_point, matched_link)
        snapped.parentRef = best_point.parentRef
        # print "For sampleID: ", probe_point.sampleID
        # print "We matched link: ", probe_point.linkPVID
        #write_points_csv(candidate_points, "points.csv")

        # link matching (updating the fields)
        # matched_link.matchedProbePoints.append(probe_point)
        control_point_index = matched_link.shapeInfo.index(best_point)
        matched_link.elevSum[control_point_index] += probe_point.point3D.elevation
        matched_link.elevCount[control_point_index] += 1

        probe_point.linkPVID = matched_link.linkPVID
        probe_point.direction = directions[best_score_index]
        #probe_point.snappedControlPoint = best_point



    # print "For sampleID: ", probe_point.sampleID
    # print "We matched link: ", probe_point.linkPVID
    #writer.writerow(probe_point.to_csv_list())
    counter = counter + 1

in_file.close()

for link in links:
    link.estimate_slope_info()
average, median, maximum, max_link = get_total_error(links)
print "Average error: ", average
print "Median error: ", median
print "Maximum error: ", maximum
print "Maximum linkID: ", max_link


#out_file.close()

# write_link_info_csv(links, "link_probe_info.csv")
# road slope derivation and evaluation
# print "Calculating Road Slope"
#
# for link in links:
#     # road slope derivation
#     if len(link.matchedProbePoints) < 2:
#         print "No point matched to this link"
#     else:
#         pass
#         # calculate road slope
#     # road slope evaluation
#     if link.slopeInfo is None:
#         print "No slope info to compare on this link"


#probe_list = assign_nodes(will_path + points_filename, will_path + links_filename)
# links = [probe.matchedLink for probe in points]

# print "Snapped points: ", len(set(best_points))
# print "Links: ", len(set(links))
# write_probe_csv(points, "probes.csv")
# write_points_csv(snapped_points, "snapped.csv")
# write_points_csv(all_candidates, "points.csv")
# write_points_csv(best_points, "matched.csv")
# write_links_csv(links, "links.csv")

# write_probe_csv(points, "probes.csv")
#write_links_csv([probe.matchedLink for probe in probe_list], "links.csv")



# parse_point_from_string('51.4965800/9.3862299/|51.4994700/9.3848799/')
