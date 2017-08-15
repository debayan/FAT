TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += core                                      #< GUI and core libraries
SUBDIRS += surfacemesh/surfacemesh                   #< a new supported model 
SUBDIRS += surfacemesh/surfacemesh_io_obj            #< loads models from file
SUBDIRS += surfacemesh/surfacemesh_render_flat       #< renders file with a style
SUBDIRS += surfacemesh/surfacemesh_filter_normalize  #< process the file with a style
