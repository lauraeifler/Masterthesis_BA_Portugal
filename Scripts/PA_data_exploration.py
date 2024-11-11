#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:08:29 2024

@author: laura
"""

import pandas as pd
import numpy as np
import geopandas as gpd

import matplotlib.pyplot as plt

#%%

path_wdpa = '/Users/laura/Masterthesis/Data/WDPA/'
path_natura = '/Users/laura/Masterthesis/Data/Natura2000/'
path_inter = '/Users/laura/Masterthesis/Data/'

#%%

wdpa = gpd.read_file(path_wdpa +  'WDPA_terrestrial_clipPRT_Aug24.shp')
wdpa = wdpa.to_crs(epsg=3763)

natura = gpd.read_file(path_natura +  'Natura2000_clipPRT.shp')
natura = natura.to_crs(epsg=3763)



#wdpa_no_overlaps = gpd.overlay(wdpa, wdpa, how='union')

#%% Dissolve WDPA to remove overlaps, keeping STATUS_YR and MANG_PLAN
wdpa_dissolved = wdpa.dissolve(by=['STATUS_YR']).reset_index()
wdpa_dissolved2 = wdpa.dissolve()

# Explode the dissolved geometries (if necessary)
wdpa_exp = wdpa_dissolved2.explode(index_parts=True).reset_index(drop=True)

#%% Dissolve Natura to remove overlaps
natura_dissolved = natura.dissolve()

# Explode the dissolved geometries (if necessary)
natura_exp = natura_dissolved.explode(index_parts=True).reset_index(drop=True)

#%% Overlay operation
gpd_inters = gpd.overlay(wdpa_exp, natura_exp, how='intersection', keep_geom_type=True)
# gpd_inters2 = gpd.overlay(wdpa_dissolved, natura_dissolved, how='intersection', keep_geom_type=True)
gpd_inters.geometry.area.sum()
# gpd_inters2.geometry.area.sum()

# Display the result with desired columns
gpd_inters_short = gpd_inters[['STATUS_YR', 'MANG_PLAN', 'geometry', 'NAME']]


gpd_union = gpd.overlay(wdpa_exp, natura_exp, how='union', keep_geom_type=True)
gpd_union['geometry'].area.sum()



## Differece between data
wdpa_only = gpd.overlay(wdpa_exp, natura_exp, how='difference')
natura_only = gpd.overlay(natura_exp, wdpa_exp, how='difference')



# area sizes
total_wdpa_only_area = wdpa_only['geometry'].area.sum()
total_natura_only_area = natura_only['geometry'].area.sum()
y datad´fraketotal_intersect_area = gpd_inters['geometry'].area.sum()
total_union_area = gpd_union['geometry'].area.sum()
total_wdpa_area = wdpa_exp['geometry'].area.sum()
total_natura_area = natura_exp['geometry'].area.sum()


# calculations
total_union_area - (total_wdpa_only_area + total_natura_only_area)

total_intersect_area/total_union_area

## proportion of areas
wdpa_non_overlap_proportion = (total_wdpa_only_area / total_wdpa_area) *100
natura_non_overlap_proportion = (total_natura_only_area / total_natura_area) * 100

100 - wdpa_non_overlap_proportion
100 - natura_non_overlap_proportion

(total_wdpa_area - total_wdpa_only_area)/total_intersect_area
(total_natura_area - total_natura_only_area)/total_intersect_area


wdpa_exp.geometry.area.sum()/gpd_union['geometry'].area.sum()
natura_exp.geometry.area.sum()/gpd_union['geometry'].area.sum()
total_natura_only_area/ total_intersect_area




wdpa_dissolved.geometry.area.sum()
wdpa_dissolved2.geometry.area.sum()


proportion_wdpa = total_intersect_area / total_wdpa_area
proportion_natura2000 = total_intersect_area / total_natura_area



wdpa = gpd.read_file(path_wdpa +  'WDPA_terrestrial_clipPRT_Aug24.shp')
wdpa_repr = wdpa.to_crs(epsg=3763)
wdpa_m = wdpa_repr.unary_union

# Convert the result back to a GeoDataFrame
wdpa_merge = gpd.GeoDataFrame(geometry=[wdpa_m], crs=wdpa_repr.crs)

wdpa_merge.geometry.area


# If the result is a MultiPolygon and you want individual polygons, you can explode it
wdpa_merge_exp = wdpa_merge.explode(index_parts=True).reset_index(drop=True)

wdpa_merge_exp.geometry.area.sum()



natura = gpd.read_file(path_natura +  'Natura2000_clipPRT.shp')
natura_repr = natura.to_crs(epsg=3763)

natura_m = natura_repr.unary_union

# Convert the result back to a GeoDataFrame
natura_merge = gpd.GeoDataFrame(geometry=[natura_m], crs=natura_repr.crs)

natura_merge.geometry.area



# If the result is a MultiPolygon and you want individual polygons, you can explode it
natura_merge_exp = natura_merge.explode(index_parts=True).reset_index(drop=True)

natura_merge_exp.geometry.area.sum()



gpd_inters = gpd.overlay(wdpa_merge_exp, natura_merge_exp, how='intersection', keep_geom_type=True)
gpd_inters['total_area_m2'] = gpd_inters['geometry'].area
gpd_inters['geometry'].area.sum()


gpd_union = gpd.overlay(wdpa_merge_exp, natura_merge_exp, how='union', keep_geom_type=True)
gpd_union['geometry'].area.sum()




## Differece between data
wdpa_only = gpd.overlay(wdpa_merge_exp, natura_merge_exp, how='difference')
natura_only = gpd.overlay(natura_merge_exp, wdpa_merge_exp, how='difference')



# area sizes
total_wdpa_only_area = wdpa_only['geometry'].area.sum()
total_natura_only_area = natura_only['geometry'].area.sum()
total_intersect_area = gpd_inters['geometry'].area.sum()
total_union_area = gpd_union['geometry'].area.sum()
total_wdpa_area = wdpa_merge_exp['geometry'].area.sum()
total_natura_area = natura_merge_exp['geometry'].area.sum()


# calculations
total_union_area - (total_wdpa_only_area + total_natura_only_area)

total_intersect_area/total_union_area

## proportion of areas
wdpa_non_overlap_proportion = (total_wdpa_only_area / total_wdpa_area) *100
natura_non_overlap_proportion = (total_natura_only_area / total_natura_area) * 100

100 - wdpa_non_overlap_proportion
100 - natura_non_overlap_proportion


