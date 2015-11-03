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
from numpy import zeros_like
from numpy import linspace
#Birch 
from sklearn.cluster import Birch
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

class Storm(object):
	"""
	Main function to generate helical symmetry
	"""
	def __init__(self, storm_file, cluster_distance, min_sample, synthetic):
		self.storm_filename = storm_file[:-4]
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
		self.sfreq = None
		self.samp = None
		self.fig, self.ax = plt.subplots(4)
		#plt.subplots_adjust(left=0.25, bottom=0.25)
		self.X = None
		self.db = None
		self.plot_handle = []
		if synthetic:
			self.__read_synthetic_data()
		else:
			self.__read_data_exp()

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

	def __read_data_exp(self):
		for line in self.storm_file.readlines()[1:]:
			csv_values = line.split(",")
			frame = float(csv_values[1])
			x = float(csv_values[2])
			y = float(csv_values[3])
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
		print "Starting DBSCAN", self.cluster_distance, self.min_sample
		#db.fit(self.all_frames_xy)
		self.X = StandardScaler().fit_transform(self.all_frames_xy)
		self.db = DBSCAN(eps=self.cluster_distance, min_samples=self.min_sample).fit(self.X)
		labels = self.db.labels_


	def save_to_pdb_dbscan(self, event):
		#Saving figure first
		filname_base = self.storm_filename+"_"+str(self.cluster_distance)+"_"+str(self.min_sample)
		self.fig.savefig(filname_base+'_full_figure.png', dpi=300)
		extent = self.ax[0].get_window_extent().transformed(self.fig.dpi_scale_trans.inverted())
		self.fig.savefig(filname_base+'_main_figure.png', bbox_inches=extent, dpi=300)
		clusters = self.db.fit_predict(self.X)
		cluster_dict = {}
		point_no = 0
		for cluster in clusters:
			cluster_dict[cluster] = []
		for cluster in clusters:
			cluster_dict[cluster].append(self.all_frames_xy[point_no])
			point_no+=1
	
		out_file = open(filname_base+".pdb","w")
		out_file_cmd = open(filname_base+".cmd","w")
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

	def update_plot(self, val):
		self.cluster_distance = round(self.sfreq.val,2)
		self.min_sample = round(self.samp.val)
		self.cluster_dbscan()
		self.show_plot()
		self.fig.canvas.draw_idle()

	def show_plot(self):
		self.ax[0].clear()
		self.ax[0].set_position([0.1,0.1,0.8,0.8])
		core_samples_mask = zeros_like(self.db.labels_, dtype=bool)
		core_samples_mask[self.db.core_sample_indices_] = True
		labels = self.db.labels_
		n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
		self.ax[0].set_title('Estimated number of clusters: %d' % n_clusters_)
		# Black removed and is used for noise instead.
		unique_labels = set(labels)
		colors = plt.cm.Spectral(linspace(0, 1, len(unique_labels)))
		for k, col in zip(unique_labels, colors):
			if k == -1:
				col = 'k'

			class_member_mask = (labels == k)
			xy = self.X[class_member_mask & core_samples_mask]
			self.ax[0].plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,markeredgecolor='k', markersize=14)

			xy = self.X[class_member_mask & ~core_samples_mask]
			self.ax[0].plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,markeredgecolor='k', markersize=6)


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
	parser.add_option("-d", "--distance", dest="cluster_distance", default=24,
					  type = 'float',
					  help="Clutering distance.  Default =  2sigma")
	parser.add_option("-m", "--min", dest="min_sample", default=0.2,
			  		type = 'int',
			  		help="Minimum sample size.  Default =  None")
	parser.add_option("-s", "--synthetic", dest="synthetic",
			  		action = 'store_true',
			  		help="Synthetic data supplied.  Default =  False")

	options, args = parser.parse_args()
	storm = Storm(options.storm_file, options.cluster_distance, options.min_sample, options.synthetic)

	clusters = storm.cluster_dbscan()
	storm.show_plot()
	axcolor = 'lightgoldenrodyellow'
	#storm.ax[1] = plt.axes([0.2, 0.0, 0.65, 0.03], axisbg=axcolor)
	#storm.ax[2]  = plt.axes([0.2, 0.05, 0.65, 0.03], axisbg=axcolor)

	storm.ax[1].set_position([0.2, 0.01, 0.65, 0.03])
	storm.ax[2].set_position([0.2, 0.05, 0.65, 0.03])

	storm.sfreq = Slider(storm.ax[1], 'Max distance to neighbor', 0.0, 1.0, valinit=0.2)
	storm.samp = Slider(storm.ax[2], 'Minimum Cluster Size', 1, 60, valinit=24, valfmt='%0.0f')

	storm.sfreq.on_changed(storm.update_plot)
	storm.samp.on_changed(storm.update_plot)

	storm.ax[3].set_position([0.9, 0.01, 0.1, 0.03])
	button = Button(storm.ax[3], 'Save', color=axcolor, hovercolor='0.975')
	button.on_clicked(storm.save_to_pdb_dbscan)
	plt.show()


