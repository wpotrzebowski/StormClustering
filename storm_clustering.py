#! /usr/bin/python
"""
STORM data analysis
"""
__author__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"
__status__ = "Prototype"

#Optparse and os import
import optparse
import os
#Numpy imports
from numpy import array
from numpy import arange
from numpy.linalg import norm
from numpy import sqrt
#BioPython import
from Bio.KDTree import KDTree
#Birch 
from sklearn.cluster import Birch
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

class PointXY(object):
    def __init__(self, x, y, z, id):
        self.coord=array([x,y,z])
	self.id = id

    def get_coord(self):
        return self.coord

    def get_id(self):
	return self.id

class NeighborSearch(object):
    """Class for neighbor searching,

    This class can be used for two related purposes:

     1. To find all atoms/residues/chains/models/structures within radius
        of a given query position.
     2. To find all atoms/residues/chains/models/structures that are within
        a fixed radius of each other.

    NeighborSearch makes use of the Bio.KDTree C++ module, so it's fast.
    """
    def __init__(self, point_list, bucket_size=10):
        """Create the object.

        Arguments:

         - atom_list - list of atoms. This list is used in the queries.
           It can contain atoms from different structures.
         - bucket_size - bucket size of KD tree. You can play around
           with this to optimize speed if you feel like it.
        """
        self.point_list=point_list
        # get the coordinates
        coord_list = [a.get_coord() for a in point_list]
        # to Nx3 array of type float
        self.coords=array(coord_list).astype("f")
        assert(bucket_size>1)
        assert(self.coords.shape[1]==3)
        self.kdt=KDTree(3, bucket_size)
        self.kdt.set_coords(self.coords)

    def _get_unique_parent_pairs(self, pair_list):
        # translate a list of (entity, entity) tuples to
        # a list of (parent entity, parent entity) tuples,
        # thereby removing duplicate (parent entity, parent entity)
        # pairs.
        # o pair_list - a list of (entity, entity) tuples
        parent_pair_list=[]
        for (e1, e2) in pair_list:
            p1=e1.get_parent()
            p2=e2.get_parent()
            if p1==p2:
                continue
            elif p1<p2:
                parent_pair_list.append((p1, p2))
            else:
                parent_pair_list.append((p2, p1))
        return uniqueify(parent_pair_list)

    # Public

    def search(self, center, radius):
        """Neighbor search.
        """
        self.kdt.search(center, radius)
        indices=self.kdt.get_indices()
        n_point_list=[]
        point_list=self.point_list
        for i in indices:
            a=point_list[i]
            n_point_list.append(a)
	return n_point_list
        


class Storm(object):
    """
    Main function to generate helical symmetry
    """
    def __init__(self, storm_file, cluster_distance, min_sample):
        self.storm_file = open(storm_file)
        self.cluster_distance = cluster_distance
	self.min_sample = min_sample
        self.minX = 10000
	self.minY = 10000
	self.frame_xy = {}
	self.all_frames_xy = []
	self.clustered_dict = {}
	self.KDPointsXY = []
	self.ns = None
	#self.__read_data_birch()
	self.__read_synthetic_data()

    def __read_data(self):
	id = 0
	for line in self.storm_file.readlines()[1:]:
	    csv_values = line.split(",")
            frame = float(csv_values[0])
            x = float(csv_values[1])
            y = float(csv_values[2])
	    z = 0
            self.KDPointsXY.append(PointXY(x,y,z,id))
            self.ns = NeighborSearch(self.KDPointsXY)
	    if x<self.minX:
		self.minX = x
	    if y<self.minY:
                self.minY = y
	    if frame not in self.frame_xy.keys():
		self.frame_xy[frame] = [(x,y,id)]
	    else:
	    	self.frame_xy[frame].append((x,y,id))
	    id+=1

    def __read_data_birch(self):
        for line in self.storm_file.readlines()[1:]:
            csv_values = line.split(",")
            frame = float(csv_values[0])
            x = float(csv_values[1])
            y = float(csv_values[2])
            z = 0
            self.all_frames_xy.append([x,y])

    def __read_synthetic_data(self):
	for line in self.storm_file.readlines():
	    csv_values = line.split(",")
            x = float(csv_values[0])
            y = float(csv_values[1])
            self.all_frames_xy.append([x,y])


    def cluster_mbk(self):
	mbk = MiniBatchKMeans(init='k-means++', n_clusters=40, batch_size=100,
                      n_init=100, max_no_improvement=10, verbose=0,
                      random_state=0)
	mbk.fit(self.all_frames_xy)
	clusters = mbk.predict(self.all_frames_xy)
	return clusters

    def cluster_birch(self):
	print "Starting Birch clustering"
	brc = Birch(branching_factor=10, n_clusters=40, threshold=self.cluster_distance,compute_labels=False)
	brc.fit(self.all_frames_xy)
	clusters = brc.predict(self.all_frames_xy)
	return clusters

    def cluster_dbscan(self):
	print "Starting DBSCAN"
	db = DBSCAN(eps=self.cluster_distance, min_samples=self.min_sample)
	#db.fit(self.all_frames_xy)
	X = StandardScaler().fit_transform(self.all_frames_xy)
	clusters = db.fit_predict(X)
	labels = db.labels_
	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	print('Estimated number of clusters: %d' % n_clusters_)
        return clusters


    def __get_distance(self, xy_cen, xy_new):
	#return norm(array(xy_cen), array(xy_new))
	return sqrt((xy_cen[0]-xy_new[0])**2+(xy_cen[1]-xy_new[1])**2)

    def cluster(self, frame_xy, cluster_distance):
	"""
	Sequential clustering
	Using neighbor search algorithm
	"""    
	counter = 1
        cluster_no = 1
	clustered_ids = []
	clustered_dict = {}
	L_r = 0.0
	cutoff_distance = 0.0
	for frame in frame_xy.keys():
	    print "Frame: ", counter
	    for xy in frame_xy[frame]:
		if not xy[2] in clustered_ids:    
		    cen_coords = array([xy[0],xy[1],0],dtype=float)
		    clustered_dict[cluster_no] = [xy]
		    for cluster_distance in arange(100.0,cluster_distance,50.0):
		    	clustered_pts = self.ns.search(cen_coords,cluster_distance)
			#Testing if it grows 
			if sqrt(len(clustered_pts)/3.14)<L_r:
			    break
			else:
			    L_r = sqrt(len(clustered_pts)/3.14)
			    cutoff_distance = cluster_distance	
		    print "Found %i at the distance %d" % (len(clustered_pts),cutoff_distance)
		    if len(clustered_pts) > 0:
		    	for clustered_pt in clustered_pts:
			    #Avoiding redudancy
			    if clustered_pt.get_id() not in clustered_ids:
			    	clustered_ids.append(clustered_pt.get_id())
				clustered_dict[cluster_no].append((clustered_pt.get_coord()[0],clustered_pt.get_coord()[1],clustered_pt.get_id()))
		    cluster_no+=1
	    counter+=1
	return clustered_dict

    def get_data(self):
	scene = display()
    	scene.fullscreen = 0
    	scene.autocenter = 1
	#scene.range = (1000,1000,10)
	self.minX-=30
	self.minY-=30
	print self.minX, self.minY
	for frame in self.frame_xy.keys():
		for xy in self.frame_xy[frame]:
		    x = xy[0]
		    y = xy[1]
		    z = 0
		    ball = sphere (pos=(x,y,z), color = color.green, radius = 30)

    def save_to_pdb_birch(self, pdb_name, clusters):
	#Sorting 
	cluster_dict = {}
	point_no = 0
	for cluster in clusters:
		cluster_dict[cluster] = []
	for cluster in clusters:
		cluster_dict[cluster].append(self.all_frames_xy[point_no])
		point_no+=1
	
	out_file = open(str(self.cluster_distance)+pdb_name,"w")
        out_file_cmd = open(str(self.cluster_distance)+".cmd","w")
        write_lines = []
        line_count = 1
        out_file_cmd.write("set bg_color white\n")
        out_file_cmd.write("vdwdefine 30\n")
        out_file_cmd.write("rep sphere\n")
	point_no = 0
	current_cluster = 0
	res_count = 1 
	for cluster in cluster_dict.keys():
	    write_lines.append("MODEL "+str(cluster+1)+"\n")
	    if cluster>=0:
		out_file_cmd.write("molmap #0."+str(cluster+1)+" 150\n")
                out_file_cmd.write("volume #0."+str(cluster+1)+" transparency 0.3\n")
		#if len(cluster_dict[cluster]) > 100:
                #	out_file_cmd.write("molmap #0."+str(cluster+1)+" 150\n")
                #	out_file_cmd.write("volume #0."+str(cluster+1)+" transparency 0.3\n")
	    res_count = 1
	    for xy_point in cluster_dict[cluster]:	
		line = "%-6s%5d  %-4s%3s %1s%4d    %8.2f%8.2f%8.2f%6.2f%6.2f\n" % ("ATOM",line_count,"C","CEN","A",res_count,xy_point[0],xy_point[1],0.00,1.0,1.0)
                write_lines.append(line)
		res_count+=1
	    write_lines.append("ENDMDL\n")
	out_file.writelines(write_lines)
        out_file.close()
        out_file_cmd.close()

 
    def save_to_pdb(self, pdb_name, clustered_dict):
	"""
	Saving clustered coordinates to fake pdb file
	"""
	out_file = open(str(self.cluster_distance)+pdb_name,"w")
	out_file_cmd = open(str(self.cluster_distance)+".cmd","w")
	write_lines = []
	line_count = 1
	out_file_cmd.write("set bg_color white\n")
	out_file_cmd.write("vdwdefine 30\n")
	out_file_cmd.write("rep sphere\n")
	out_file_cmd.write("col red @CA\n")
	for centroid in clustered_dict.keys():
            res_count = 1
	    write_lines.append("MODEL "+str(centroid)+"\n")
	    real_centroid = array([0,0],dtype=float)
            for clustered in clustered_dict[centroid]:
            	real_centroid +=array([clustered[0],clustered[1]],dtype=float)
            	line = "%-6s%5d  %-4s%3s %1s%4d    %8.2f%8.2f%8.2f%6.2f%6.2f\n" % ("ATOM",line_count,"C","CEN","A",res_count,clustered[0],clustered[1],0.00,1.0,1.0)
		write_lines.append(line)
            real_centroid /=len(clustered_dict[centroid])
            line = "%-6s%5d  %-4s%3s %1s%4d    %8.2f%8.2f%8.2f%6.2f%6.2f\n" % ("ATOM",line_count,"CA","CEN","A",res_count, real_centroid[0],real_centroid[1],0.00,1.0,1.0)
            write_lines.append(line)
        write_lines.append("ENDMDL\n")
        out_file_cmd.write("molmap #0."+str(centroid)+" 150\n")
        out_file_cmd.write("volume #0."+str(centroid)+" transparency 0.3\n")
        line_count+=1
	out_file.writelines(write_lines)
	out_file.close()
	out_file_cmd.close()

if '__main__' == __name__:
    doc = """
    Clusters storm data in a sequntial manner
    usage: python storm_cluster.py --help
    """
    print doc
    usage = "usage: %prog [options] args"
    option_parser_class = optparse.OptionParser
    parser = option_parser_class( usage = usage, version='0.1' )

    parser.add_option("-f", "--file", dest="storm_file", default = None,
                      type = 'string',
                      help="STORM file in a csv format from ThunderSTORM [OBLIGATORY]")
    parser.add_option("-d", "--distance", dest="cluster_distance", default =None,
                      type = 'float',
                      help="Clutering distance.  Default =  2sigma")
    parser.add_option("-m", "--min", dest="min_sample", default =None,
		      type = 'int',
		      help="Minimum sample size.  Default =  None")


    options, args = parser.parse_args()
    storm = Storm(options.storm_file, options.cluster_distance, options.min_sample)
    clusters = storm.cluster_dbscan()
    print clusters
    #First iteration
    #cluster_dict = storm.cluster(storm.frame_xy, options.cluster_distance)
    storm.save_to_pdb_birch( "output.pdb", clusters)
        

