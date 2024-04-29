######CORAL DIVERSITY
####RASTERIZE MARINE BIOGEOREGIONS
#from shapely.wkb import loads
import rasterio
import matplotlib.pyplot as plt
import csv
import os
import statsmodels
import statsmodels.api as sm
import pandas as pd
from numpy import array,where,zeros,linspace,meshgrid,nan,isnan,ma,corrcoef,mean
from osgeo import gdal,gdal_array
from collections import Counter
from itertools import combinations
from mpl_toolkits.basemap import Basemap
from sklearn.ensemble import RandomForestRegressor



def pixel2coord(x, y, xoff, a, b, yoff, d, e):
	xp = a * x + b * y + a * 0.5 + b * 0.5 + xoff
	yp = d * x + e * y + d * 0.5 + e * 0.5 + yoff
	return(xp, yp)


def point_to_cell(point_x, point_y, cellx, celly, xmin, ymax):
    col = int((point_x - xmin) / cellx)
    row = int((point_y - ymax) / -celly)
    return row,col



susc = list(csv.reader(open('susc.csv','r')))
len(set([i[0] for i in susc[1:] if i[0]!='']))
len(set([i[1] for i in susc[1:] if i[0]!='']))
len(set([i[2] for i in susc[1:] if i[0]!='']))

zoo_f = list(csv.reader(open('GeoSymbio.csv','r')))
 for row in zoo_f:
  	if row[10] == 'Montastraea':
  		row[9] = 'Montastraeidae'


ok_species = set([i[0] for i in susc[1:]]) & set([i[12] for i in zoo_f])


data = []
for sp in ok_species:
	row = []
	for i in susc:
		if i[0]==sp:
			row.append(list(map(int,i[3:6]))+i[1:3])
	row.sort()
	prev = row[-1][1]/row[-1][0]
	mort = row[-1][2]/row[-1][0]
	gen,fam = row[0][3:5]
	data.append([sp,gen,fam,prev,mort])





clades = sorted(list(set([i[0] for i in zoo_f[1:]])))

#genera in geosymbio

susc_sp = [[i[0],i[3],i[4],1,0,0] for i in data]
all_gen = set([i[1] for i in data])
all_fam = set([i[2] for i in data])
susc_gen = []
for i in all_gen:
	row_prev,row_mort = [],[]
	for j in data:
		if j[1]==i:
			row_prev.append(j[3])
			row_mort.append(j[4])
	susc_gen.append([i,mean(row_prev),mean(row_mort),0,1,0])



susc_fam = []
for i in all_fam:
	row_prev, row_mort = [], []
	for j in data:
		if j[2] == i:
			row_prev.append(j[3])
			row_mort.append(j[4])
	susc_fam.append([i, mean(row_prev), mean(row_mort),0,0,1])


susc = susc_sp+susc_gen+susc_fam

susc_ = []
zoo_susc = []
for i in range(len(susc)):
	row = []
	for j in zoo_f:
		if susc[i][0] in j[9:13]:
			row.append(j[0])
	if row!=[]:
		susc_.append(susc[i])
		zoo_susc.append(row)



len([i for i in susc_ if i[-3]==1])
len([i for i in susc_ if i[-2]==1])
len([i for i in susc_ if i[-1]==1])

zoo_susc = [[Counter(i)[j]/float(len(i)) for j in clades] for i in zoo_susc]

x = [zoo_susc[i]+susc_[i][-3:] for i in range(len(susc_))]
y_prev = [float(i[1]) for i in susc_]
y_mort = [float(i[2]) for i in susc_]

out = open('susc_table.csv','w')
out.write(','.join(['taxon','prevalence','mortality','sp','gen','fam','a','b','c','d','e','f','g','h','i'])+'\n')
for i in range(len(susc_)):
	out.write(','.join(map(str,susc_[i]+zoo_susc[i]))+'\n')


out.close()

df_prev,df_mort = [],[]
for i in range(len(x)):
	df_prev.append([i,y_prev[i]]+x[i])
	df_mort.append([i,y_mort[i]]+x[i])


df_prev = pd.DataFrame(df_prev,columns=['id','prevalence','a','b','c','d','e','f','g','h','i','sp','gen','fam'])
df_mort = pd.DataFrame(df_mort,columns=['id','mortality','a','b','c','d','e','f','g','h','i','sp','gen','fam'])

#test out of bag prev
cols = ['a','b','c','d','e','f','g','h','i']
best_mod_prev = RandomForestRegressor(n_estimators=1000, oob_score=False)
obs,pred = [],[]
for lo in range(len(x)):
	df_train = df_prev.loc[df_prev['id']!=lo]
	df_test = df_prev.loc[df_prev['id']==lo]
	X = pd.DataFrame(df_train,columns=list(cols)+['sp','gen','fam'])
	y = pd.DataFrame(df_train,columns=['prevalence'])
	X_test = pd.DataFrame(df_test,columns=list(cols)+['sp','gen','fam'])
	y_test = pd.DataFrame(df_test,columns=['prevalence'])
	best_mod_prev.fit(X, y.values.ravel())
	pred.append([best_mod_prev.predict(X_test)[0]]+x[lo][-3:])
	obs.append([y_test.values[0][0]]+x[lo][-3:])


r2_prev = corrcoef([i[0] for i in obs],[i[0] for i in pred])[0][1]**2
r2_sp_prev = corrcoef([i[0] for i in obs if i[1]==1],[i[0] for i in pred if i[1]==1])[0][1]**2
r2_gen_prev = corrcoef([i[0] for i in obs if i[2]==1],[i[0] for i in pred if i[2]==1])[0][1]**2
r2_fam_prev = corrcoef([i[0] for i in obs if i[3]==1],[i[0] for i in pred if i[3]==1])[0][1]**2




X = pd.DataFrame(df_prev,columns=list(cols)+['sp','gen','fam'])
y = pd.DataFrame(df_prev,columns=['prevalence'])
mod_prev = RandomForestRegressor(n_estimators=1000, oob_score=True)
mod_prev.fit(X, y.values.ravel())
mod_prev.oob_score_

###mortality
best_mod_mort = RandomForestRegressor(n_estimators=1000, oob_score=False)
obs,pred = [],[]
for lo in range(len(x)):
	df_train = df_mort.loc[df_mort['id']!=lo]
	df_test = df_mort.loc[df_mort['id']==lo]
	X = pd.DataFrame(df_train,columns=list(cols)+['sp','gen','fam'])
	y = pd.DataFrame(df_train,columns=['mortality'])
	X_test = pd.DataFrame(df_test,columns=list(cols)+['sp','gen','fam'])
	y_test = pd.DataFrame(df_test,columns=['mortality'])
	best_mod_mort.fit(X, y.values.ravel())
	pred.append([best_mod_mort.predict(X_test)[0]]+x[lo][-3:])
	obs.append([y_test.values[0][0]]+x[lo][-3:])



r2_mort = corrcoef([i[0] for i in obs],[i[0] for i in pred])[0][1]**2
r2_sp_mort = corrcoef([i[0] for i in obs if i[1]==1],[i[0] for i in pred if i[1]==1])[0][1]**2
r2_gen_mort = corrcoef([i[0] for i in obs if i[2]==1],[i[0] for i in pred if i[2]==1])[0][1]**2
r2_fam_mort = corrcoef([i[0] for i in obs if i[3]==1],[i[0] for i in pred if i[3]==1])[0][1]**2


X = pd.DataFrame(df_mort,columns=list(cols)+['sp','gen','fam'])
y = pd.DataFrame(df_mort,columns=['mortality'])

mod_mort = RandomForestRegressor(n_estimators=1000, oob_score=True)
mod_mort.fit(X, y.values.ravel())
mod_mort.oob_score_

####fit

print (round(r2_prev,2))
print (round(r2_sp_prev,2))
print (round(r2_gen_prev,2))
print (round(r2_fam_prev,2))


print (round(r2_mort,2))
print (round(r2_sp_mort,2))
print (round(r2_gen_mort,2))
print (round(r2_fam_mort,2))



####
zoo_a = [i for i in zoo_f if i[7]=='Anthozoa']
rst_fn = 'reef.tif'
rst = rasterio.open(rst_fn)
meta = rst.meta.copy()
meta.update(compress='lzw')
reef_mat = gdal_array.LoadFile(rst_fn)
sss=reef_mat.shape
occs = where(reef_mat>0)
reef_cells = zip(occs[0],occs[1])


div_mat = zeros(sss)
meas_1, meas_3 = [zeros(sss) for i in range(2)]
gen = sorted(list(set([i[10].replace(' ','') for i in zoo_a])))
net_gen = [[i,[]] for i in gen]
for i in zoo_a:
	net_gen[gen.index(i[10].replace(' ',''))][1].append(i[0])


net_gen_c = [[i[0],[Counter(i[1])[j]/float(len(i[1])) for j in clades]+[0,1,0]] for i in net_gen]
susc_dict_gen = []
for i in net_gen_c:
	pred_prev = mod_prev.predict(pd.DataFrame([i[1]],columns=cols+['sp','gen','fam']))[0]
	pred_mort = mod_mort.predict(pd.DataFrame([i[1]],columns=cols+['sp','gen','fam']))[0]
	susc_dict_gen.append([i[0],[pred_prev,pred_mort]])


susc_dict_gen = dict(susc_dict_gen)

#############################
###Derive abundance vs. diversity risk
caribbean = list(csv.reader(open('caribbean_transects.csv','r')))
caribbean_genera = caribbean[0][7:]

caribbean_genera_prev_mort = [susc_dict_gen.get(i,[0,0]) for i in caribbean_genera]
caribbean_genera_prev = array([i[0] for i in caribbean_genera_prev_mort])
caribbean_genera_mort = array([i[1] for i in caribbean_genera_prev_mort])
caribbean_genera_prev_per_mort = caribbean_genera_prev*caribbean_genera_mort

caribbean_res = []
for i in caribbean[1:]:
	genera_cover = array(list(map(float, i[7:])))
	genera_pa = genera_cover > 0
	genera_n = genera_pa.sum()
	genera_cover_frac = genera_cover / (genera_cover).sum()
	prev_pres_abs = genera_pa * caribbean_genera_prev
	prev_per_mort_pres_abs = genera_pa * caribbean_genera_prev_per_mort
	prev_cover = genera_cover_frac * caribbean_genera_prev
	prev_per_mort_cover = genera_cover_frac * caribbean_genera_prev_per_mort
	###new_scores
	wI = (caribbean_genera_prev * genera_n * genera_cover_frac).sum()
	wIII = (caribbean_genera_prev_per_mort * genera_n * genera_cover_frac).sum()
	I = prev_pres_abs.sum()
	II = I / genera_n
	III = prev_per_mort_pres_abs.sum()
	IV = III / genera_n
	### weighted scores
	I_c = prev_cover.sum()
	II_c = I_c / genera_n
	III_c = prev_per_mort_cover.sum()
	IV_c = III_c / genera_n
	caribbean_res.append([I, I_c, II, II_c, III, III_c, IV, IV_c, wI, wIII, genera_n, genera_cover.sum(), genera_cover_frac.std()])
	print(i, len(caribbean_res))




out = open('caribbean_metrics.csv','w')
out.write('I,I_c,II,II_c,III,III_c,IV,IV_c,wI,wIII,coral_genera_n,total_coral_cover,coral_cover_sd\n')
for i in caribbean_res:
    out.write(','.join(map(str,i))+'\n')

out.close()


prev_mats,imp_mats = [],[]
for i in susc_dict_gen.keys():
	prev,mort = susc_dict_gen[i]
	try:
		mat = gdal_array.LoadFile('./coral_ranges_reef/'+i+'.tif')
		div_mat += mat
		if prev!='na':
			prev_mats.append(mat*prev)
			imp_mats.append(mat*prev*mort)
	except:
		print (i,prev)


for m in range(len(prev_mats)):
	meas_1 += prev_mats[m]
	meas_3 += imp_mats[m]


meas_2 = meas_1/div_mat
meas_4 = meas_3/div_mat


meas_2*= 0.8629 #mean prevalence
meas_1*=0.7826 #cumulative prevalence
meas_4*=0.7202 #mean mortality
meas_3*=0.5142 #cumulative mortality




cm = plt.get_cmap('plasma')

names = ['cumulative_prevalence','mean_prevalence','cumulative_mortality','mean_mortality']
sc = 0
for data in  [meas_1,meas_2,meas_3,meas_4]:
	data[where(data==0)] = nan
	plt.rcParams['hatch.linewidth'] = 0.001
	m = Basemap(resolution='l', projection='cyl',llcrnrlat=-45, llcrnrlon=15, urcrnrlat=45, urcrnrlon= 310)
	m.drawcoastlines(linewidth=0.1)
	m.fillcontinents(color='darkgrey', lake_color='darkgrey')#,zorder=0)  # zorder=0 to paint over continents
	data_m = ma.masked_where(isnan(data),data)
	x = linspace(-180, 180, data.shape[1])
	y = linspace(35, -35, data.shape[0])
	xx, yy = meshgrid(x, y)
	xx[where(xx<0)]+=360
	im = m.pcolormesh(xx, yy, data_m, cmap=cm)
	plt.colorbar(orientation = 'horizontal',pad=0.01)
	plt.savefig(names[sc]+'_map.png', dpi=300,bbox_inches='tight')
	clear = [plt.clf() for i in range(10000)]
	sc += 1



