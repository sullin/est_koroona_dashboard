#!/usr/bin/env python3
# Koroona dashboard InfluxDB ja Grafana'ga.
# MIT License
# Copyright (c) 2021 Anti Sullin

import urllib.request
import json
import locale
import socket
import datetime
import numpy, scipy
import math
from scipy.optimize import curve_fit
from scipy.stats import norm
from influxdb import InfluxDBClient

# Influxdb server
influxhost = '192.168.32.2'
influxport = 8086

# Regioonid ja nende elanike arvud
regions = {
	"Harju":	598059,
	"Hiiu":	9387,
	"Ida-Viru":	136240,
	"Järva":	28734,
	"Jõgeva":	30286,
	"Lääne":	20507,
	"Lääne-Viru":	59325,
	"Pärnu":	25006,
	"Põlva":	85938,
	"Rapla":	33311,
	"Saare":	33108,
	"Tartu":	152977,
	"Valga":	28370,
	"Viljandi":	46371,
	"Võru":	35782,
	"Tundmatu": 0,
	"Välismaa": 0,
	"Kokku":	1323401
}


# Arvutab viimase "days" päeva juurdekasvu
def calc_incr(total, days):
	ret = {}
	for d in total:
		dref = d - datetime.timedelta(days=days)
		if not dref in total:
			continue
		ret[d]={}
		for mk in total[d]:
			ret[d][mk]=total[d][mk]-total[dref][mk]
	return ret

# Arvutab igas piirkonnas normaliseeritud arvu "cnt" elaniku kohta
def calc_norm(total, regions, cnt):
	ret = {}
	for d in total:
		ret[d]={}
		for mk in total[d]:
			if not mk in regions or regions[mk]==0:
				continue
			ret[d][mk] = cnt * total[d][mk] // regions[mk]
	return ret


# Arvutab n-päeva kordajat (R0), kasutades logaritmilist regressiooni; vaatab "window" pikkust akent
def calc_logest(pos, window, transmission_delay):
	# https://stackoverflow.com/questions/55780106/how-to-mimic-excels-logest-function-in-python
	# LOGEST from https://support.office.com/en-us/article/logest-function-f27462d8-3657-4030-866b-a272c1d18b4b
	def func(x, b, m):
		y = b * m**x
		return y
	x=numpy.array(range(window))
	ret={}
	
	for d in pos:
		ret[d]={}
		for mk in pos[d]:
			try:
				y=[]
				for i in range(window):
					dd = d-datetime.timedelta(days=window-1-i)
					v=0
					if dd in pos and mk in pos[dd]:
						v=pos[dd][mk]
					if v == 0:
						v = 0.0001
					y.append(v)
				
				# https://stackoverflow.com/questions/3433486/how-to-do-exponential-and-logarithmic-curve-fitting-in-python-i-found-only-poly
				y = numpy.array(y)
				fittedParameters = numpy.polyfit(x, numpy.log(y), 1, w=numpy.sqrt(y))
				val = (numpy.e**fittedParameters[0])**transmission_delay
				if val > 1000 or val < 1/1000:
					val = 1
				ret[d][mk]=val
			except:
				pass
	return ret

# Nakkuse saamine on ca 1 nädal enne positiivset diagnoosi.
# Graafik lõppeb offs+w päeva minevikus.
# https://www.youtube.com/watch?v=LnQcbAKWkPE&t=487s
def calc_exposure(daily, offs, w):
	coefs=[]
	for i in range(2*w+1):
		coefs.append(norm.pdf(i-w,scale=3.0))

	start = min(daily.keys())+datetime.timedelta(days=w)

	ret = {}
	for d in daily:
		if d < start:
			continue
		dout = d - datetime.timedelta(days=offs+w)
		ret[dout] = {}
		for mk in daily[d]:
			v=0
			for i in range(len(coefs)):
				dd = d - datetime.timedelta(days=i)
				if not dd in daily:
					continue
				v += daily[dd][mk] * coefs[i]
			ret[dout][mk] = v
	return ret
		
# kasutame online andmeid või testime kohaliku failiga
if True:
	url="https://opendata.digilugu.ee/opendata_covid19_test_county_all.json"
	with urllib.request.urlopen(url) as f:
		data = json.loads(f.read().decode())
else:
	f = open("opendata_covid19_test_county_all.json", "r")
	data = json.loads(f.read())
	f = 0

# kokku
total={}
maakonnad=set()

# loeme sisendi, konverteerime maakonna nimed, arvutame summa
for d in data:
	if d['ResultValue'] != "P":
		continue
	kpv=datetime.datetime.strptime(d['StatisticsDate'],"%Y-%m-%d")
	if not kpv in total:
		total[kpv]={}
		total[kpv]["Kokku"]=0
	mk=d['Country']+" "+d['County']
	mk=mk.replace(" maakond","").replace("Eesti ", "").strip()
	val=d['TotalCases']
	total[kpv][mk]=val
	total[kpv]["Kokku"]+=val


# viimane ööpäev
daily=calc_incr(total,1)

# viimane ööpäev 100k kohta
daily100k=calc_norm(daily, regions, 100000)

# viimane ndl
weekly=calc_incr(total,7)

# viimane ndl 100k kohta
weekly100k=calc_norm(weekly, regions, 100000)

# hetkel haigete arv
pos=calc_incr(total,14)

# hetkel haigete arv 100k kohta
pos100k=calc_norm(pos, regions, 100000)

# r0: seda saaks ennustada erinevate andmete pealt, tuleb teha kompromiss 
# tulemuse mürasuse ja muutumiskiiruse vahelt.
# Siin arvutan regressioni viimase 14 päeva haigete arvu viimase 7 punkti pealt.
r0 = calc_logest(pos, 7, 5)
dailyinc = calc_logest(pos, 7, 1)

# haigestumise hetk
exposure = calc_exposure(daily, 7, 5)
exposure100k=calc_norm(exposure, regions, 100000)

# Väljund CSV faili
def print_csv(pr, regions):
	out=""
	for mk in regions:
		out += ","+mk
	print(out)

	for d in pr:
		out = str(d.timestamp())
		for mk in regions:
			out+=","
			if mk in pr[d]:
				out += str(pr[d][mk])
		print(out)

# Väljund InfluxDBsse
ifxclient = InfluxDBClient(host=influxhost, port=influxport)
def send_influx(data, key, tags):
	global ifxclient
	upl=[]
	for d in data:
		for mk in data[d]:
			line = "mk=%s %s=%f %d"%(mk, key, data[d][mk], int(d.timestamp())*1000000000)
			if tags:
				line = tags+","+line
			line="koroona,"+line
			upl.append(line)
	ifxclient.write_points(upl, database='grafana', time_precision=None, batch_size=10000, protocol='line')


# debug: kirjutame csv faili?
if False:
	print_csv(r0, regions)
else:
	send_influx(total, "total", None)
	send_influx(daily, "daily", None)
	send_influx(daily100k, "daily100k", None)
	send_influx(weekly, "weekly", None)
	send_influx(weekly100k, "weekly100k", None)
	send_influx(pos, "pos", None)
	send_influx(pos100k, "pos100k", None)
	send_influx(r0, "r0", None)
	send_influx(dailyinc, "dailyinc", None)
	send_influx(exposure, "exposure", None)
	send_influx(exposure100k, "exposure100k", None)


