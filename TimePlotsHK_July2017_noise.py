#!/usr/bin/env python

import copy
import os, sys, ROOT, math
from math import exp, atan, sqrt, log, atan2, tan
from array import array
from ROOT import gStyle
#from rootpy.interactive import wait

rndm = ROOT.TRandom(1)
SpeedOfLight = 30

def readLorentzVector():
  print sys.argv[1]
  f = ROOT.TFile(sys.argv[1])
  f2=ROOT.TFile(sys.argv[2])
  
  mychain2 =f2.Get("HSCPTree/tree;1") #Muons
  mychain = f.Get("HSCPTree/tree;1") #HSCP
  entries = mychain.GetEntriesFast()
  
  #Save file parameters
  str_save="plot_noise/"
  #str_save="plots/new/Go"
  #str_save="plots/HSCP_old_NoPU/"
  #str_save="plots/910pre3/"
  str_form=".pdf"
  
  c1=ROOT.TCanvas()
  seriousHistL=[]
  cseriousHistL=[]
  sRange=[[-10,100],[-10,100],[-10,100],[-10,100],[-35,85],[-10,200],[-10,200],[-0.011,0.011],[-0.2,0.2], [0,10], [8, 13], [8, 13]]*3
  numBin=150
  for i in range(12):
    s="seriousHist"+str(i)
    sc="cseriousHist"+str(i)
    seriousHistL.append(ROOT.TH1F(s,"placeholder",numBin,sRange[i][0],sRange[i][1]))
    cseriousHistL.append(ROOT.TH1F(sc,"placeholder",numBin,sRange[i][0],sRange[i][1]))
  numBin=150
  RangeSplit=[i*1000/150 for i in range(0,150)]


  def findClosestDigi (x, y, z, pos1, pos2) :
    return (x-pos1[0])**2+(y-pos1[1])**2+(z-pos1[2])**2<(x-pos2[0])**2+(y-pos2[1])**2+(z-pos2[2])**2
 

  def isWithinTheCone (x, y, z, pos1) :
    return sqrt((x-pos1[0])**2+(y-pos1[1])**2+(z-pos1[2])**2) < 10
           
          
  cMax=[-100000]*10
  cMin=[100000]*10
  
  n1_isIRPC = 0
  n2_isIRPC = 0
  
  Max=[-100000]*10
  Min=[100000]*10
  global HSCPmass
  HSCPmass = 0
  def timePlots (chain,HistL,Min,Max):



    count = -1
    for event in chain :
      if count > 1000:
        break

      count = count+1

      n=len(event.simDigi1_isIRPC)
      n2=len(event.simDigi2_isIRPC)
      m=len(event.rpcHit_isIRPC)
      pos11=(0,0,0)
      pos12=(0,0,0)
      pos21=(0,0,0)
      pos22=(0,0,0)

      r011,r012,r021,r022=[0,0,0,0]
      
      for i in range(n):
        if event.simDigi1_isIRPC[i]:
          if (abs(event.simDigi1_z[i])<1000 and abs(event.simDigi1_z[i])>900):
            pos11=(event.simDigi1_x[i],event.simDigi1_y[i],event.simDigi1_z[i])
            r011 = (event.simDigi1_tof[i]- event.simDigi1_t0[i])*SpeedOfLight
          elif abs(event.simDigi1_z[i])>1000 :
            pos12=(event.simDigi1_x[i],event.simDigi1_y[i],event.simDigi1_z[i])
            r012 = (event.simDigi1_tof[i]- event.simDigi1_t0[i])*SpeedOfLight
      for i in range(n2):
        if event.simDigi2_isIRPC[i]:
          if (abs(event.simDigi2_z[i])<1000 and abs(event.simDigi2_z[i])>900):
            pos21=(event.simDigi2_x[i],event.simDigi2_y[i],event.simDigi2_z[i])
            r021 = (event.simDigi2_tof[i]- event.simDigi2_t0[i])*SpeedOfLight
          elif abs(event.simDigi2_z[i])>1000 :
            pos22=(event.simDigi2_x[i],event.simDigi2_y[i],event.simDigi2_z[i])
            r022 = (event.simDigi2_tof[i]- event.simDigi2_t0[i])*SpeedOfLight

      t01,t02,tr1,tr2,dt,chi=[0,0,0,0,0,0]



      fakeRE31_plus = []
      fakeRE41_plus = []
      fakeRE31_minus = []
      fakeRE41_minus = []


      for i in range(m):
        if event.rpcHit_isIRPC[i]:
          if abs(event.rpcHit_z[i])<1000 and abs(event.rpcHit_z[i])>900:
             if not isWithinTheCone(event.rpcHit_x[i],event.rpcHit_y[i],event.rpcHit_z[i],pos11) and not isWithinTheCone(event.rpcHit_x[i],event.rpcHit_y[i],event.rpcHit_z[i],pos21):
               tr1=event.rpcHit_time[i]-(math.sqrt(event.rpcHit_y[i]**2+event.rpcHit_x[i]**2+event.rpcHit_z[i]**2)-r011)/SpeedOfLight
               if event.rpcHit_z[i] > 0: 
                 fakeRE31_plus.append(i)
               else : 
                 fakeRE31_minus.append(i)


          elif abs(event.rpcHit_z[i])>1000:
            if not isWithinTheCone(event.rpcHit_x[i],event.rpcHit_y[i],event.rpcHit_z[i],pos21) and not isWithinTheCone(event.rpcHit_x[i],event.rpcHit_y[i],event.rpcHit_z[i],pos22):


#######################
              if event.rpcHit_z[i] > 0: 
                fakeRE41_plus.append(i)
              else : 
                fakeRE41_minus.append(i)

      l31_plus=len(fakeRE31_plus)
      l31_minus=len(fakeRE31_minus)

      l41_plus=len(fakeRE41_plus)
      l41_minus=len(fakeRE41_minus)

      for i in range(l31_plus):
        for j in range(l41_plus):
          ii = fakeRE31_plus[i]
          jj = fakeRE41_plus[j]
          tr1=event.rpcHit_time[ii]-(math.sqrt(event.rpcHit_y[ii]**2+event.rpcHit_x[ii]**2+event.rpcHit_z[ii]**2)-r012)/SpeedOfLight
          tr2=event.rpcHit_time[jj]-(math.sqrt(event.rpcHit_y[jj]**2+event.rpcHit_x[jj]**2+event.rpcHit_z[jj]**2)-r012)/SpeedOfLight
          
          xi = event.rpcHit_x[ii]
          xj = event.rpcHit_x[jj]

          yi = event.rpcHit_y[ii]
          yj = event.rpcHit_y[jj]        

          zi = event.rpcHit_z[ii]
          zj = event.rpcHit_z[jj]   

          ri = sqrt(xi*xi+yi*yi)
          rj = sqrt(xj*xj+yj*yj)

          phii = atan2(yi, xi)
          phij = atan2(yj, xj)

          etai = - log(tan(atan2(ri, zi)/2))
          etaj = - log(tan(atan2(rj, zj)/2))

      #    print etai, " ", etaj

          if abs(phij-phii) < 0.1 and abs(etai - etaj) < 0.1:

            dt=tr2-tr1
            HistL[4].Fill(dt)
            
          



#      print "size31 plus ",len(fakeRE31_plus)
#      print "size41 plus  ",len(fakeRE41_plus)

#      print "size31 minus ",len(fakeRE31_minus)
#      print "size41 minus  ",len(fakeRE41_minus)


  timePlots(mychain,seriousHistL,Min,Max)
  timePlots(mychain2,cseriousHistL,Min,Max)


  c1=ROOT.TCanvas()
  cseriousHistL[4].SetLineColor(ROOT.kRed)
  cseriousHistL[4].GetXaxis().SetTitle('Time delay difference')
  cseriousHistL[4].GetXaxis().SetRangeUser(-50, 100)
  cseriousHistL[4].GetYaxis().SetTitle('N_Particles')
  cseriousHistL[4].SetTitle('Difference between second and first chamber')
  seriousHistL[4].SetLineColor(ROOT.kBlue)

  norm_mu = cseriousHistL[4].Integral()  # MG normalize to unity
  cseriousHistL[4].Scale(1./norm_mu)  # MG normalize to unity
  norm_hscp = seriousHistL[4].Integral()  # MG normalize to unity
  seriousHistL[4].Scale(1./norm_hscp)  # MG normalize to unity


  cseriousHistL[4].Draw("HIST")
  seriousHistL[4].Draw("SAME HIST")
  legend = ROOT.TLegend(0.13,0.70,0.30,0.87)
  legend.AddEntry(cseriousHistL[4],"Muons","l")
  legend.AddEntry(seriousHistL[4],"HSCP","l")
  legend.Draw() 

  cseriousHistL[4].SetMaximum(0.3)

  c1.SaveAs(str_save+"timeDiff"+str_form)



readLorentzVector()

