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
from numpy import zeros_like
from numpy import linspace
#BioPython import
from Bio.KDTree import KDTree
#Birch 
from sklearn.cluster import Birch
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

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
		#db.fit(self.all_frames_xy)
		X = StandardScaler().fit_transform(self.all_frames_xy)
		db = DBSCAN(eps=self.cluster_distance, min_samples=self.min_sample).fit(X)
		self.__show_plot(db,X)
		clusters = db.fit_predict(X)

		return clusters


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

		out_file_cmd.write("set bg_color white\n")
		out_file_cmd.write("vdwdefine 30\n")
		out_file_cmd.write("rep sphere\n")

		for cluster in cluster_dict.keys():
			write_lines.append("MODEL "+str(cluster+1)+"\n")
			if cluster>=0:
				out_file_cmd.write("molmap #0."+str(cluster+1)+" 150\n")
				out_file_cmd.write("volume #0."+str(cluster+1)+" transparency 0.3\n")
			res_count = 1
			line_count = 1
			for xy_point in cluster_dict[cluster]:
				line = "%-6s%5d  %-4s%3s %1s%4d    %8.2f%8.2f%8.2f%6.2f%6.2f\n" % ("ATOM",line_count,"C","CEN","A",res_count,xy_point[0],xy_point[1],0.00,1.0,1.0)
				write_lines.append(line)
				res_count+=1
				line_count+=1
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

	def __show_plot(self, db, X):
		core_samples_mask = zeros_like(db.labels_, dtype=bool)
		core_samples_mask[db.core_sample_indices_] = True
		labels = db.labels_
		n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
		# Black removed and is used for noise instead.
		unique_labels = set(labels)
		colors = plt.cm.Spectral(linspace(0, 1, len(unique_labels)))
		for k, col in zip(unique_labels, colors):
			if k == -1:
				col = 'k'

			class_member_mask = (labels == k)
			xy = X[class_member_mask & core_samples_mask]
			plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,markeredgecolor='k', markersize=14)

			xy = X[class_member_mask & ~core_samples_mask]
			plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,markeredgecolor='k', markersize=6)

		plt.title('Estimated number of clusters: %d' % n_clusters_)

		axcolor = 'lightgoldenrodyellow'
		axfreq = plt.axes([0.2, 0.0, 0.65, 0.03], axisbg=axcolor)
		axamp  = plt.axes([0.2, 0.05, 0.65, 0.03], axisbg=axcolor)

		a0 = 24
		f0 = 0.5
		sfreq = Slider(axfreq, 'Max distance to neighbor', 0.0, 1.0, valinit=f0)
		samp = Slider(axamp, 'Minimum Cluster Size', 1, 60, valinit=a0)

		plt.show()

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


