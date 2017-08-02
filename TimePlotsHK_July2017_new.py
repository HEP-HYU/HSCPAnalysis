#!/usr/bin/env python

import copy
import os, sys, ROOT, math
from math import exp, atan, sqrt, log, atan2, tan
from array import array
from ROOT import gStyle
#from rootpy.interactive import wait

rndm = ROOT.TRandom(1)
SpeedOfLight = 30

def readLorentzVector(Pileup):
  print sys.argv[1]
  f = ROOT.TFile(sys.argv[1])
  f2=ROOT.TFile(sys.argv[2])
  
  mychain2 =f2.Get("HSCPTree/tree;1") #Muons
  mychain = f.Get("HSCPTree/tree;1") #HSCP
  entries = mychain.GetEntriesFast()
  
  #Save file parameter
  #str_save="plots/new/release/001/1599/"
  #str_save = "tmp/"
  #str_save="plot_HSCPLGW25_noise_BgPU0/"
  str_save = "plot_NOnoise_noPU/"
  #str_save="plots/new/Go"
  #str_save="plots/HSCP_old_NoPU/"
  #str_save="plots/910pre3/"
  #str_form=".root"
  str_form=".pdf"
  
  c1=ROOT.TCanvas()
  seriousHistL=[]
  cseriousHistL=[]
  sRange=[[-10,100],[-10,100],[-10,100],[-10,100],[-35,85],[-10,200],[-10,200],[-0.011,0.011],[-0.2,0.2], [0,10], [8, 13], [8, 13], [0, 3], [0, 3]]*3
  numBin=[     50,      50,        50,      50,        50,      50,      50,          100,          100,       50,     50,      50,     100,      100]
  for i in range(14):
    s="seriousHist"+str(i)
    sc="cseriousHist"+str(i)
    seriousHistL.append(ROOT.TH1F(s,"placeholder",numBin[i],sRange[i][0],sRange[i][1]))
    cseriousHistL.append(ROOT.TH1F(sc,"placeholder",numBin[i],sRange[i][0],sRange[i][1]))

  s="seriousHist14"
  sc="cseriousHist14"
  seriousHistL.append(ROOT.TH2F(s,"placeholder",200,900,1100, 200, 200, 400))
  cseriousHistL.append(ROOT.TH2F(sc,"placeholder",200,900,1100, 200, 200, 400))

  s="seriousHist15"
  sc="cseriousHist15"
  seriousHistL.append(ROOT.TH2F(s,"placeholder",200, 200, 400, 200, 200, 400))
  cseriousHistL.append(ROOT.TH2F(sc,"placeholder",200, 200, 400, 200, 200, 400))
  


  def findClosestDigi (x, y, z, pos1, pos2) :
    return (x-pos1[0])**2+(y-pos1[1])**2+(z-pos1[2])**2<(x-pos2[0])**2+(y-pos2[1])**2+(z-pos2[2])**2


  def isWithinTheCone (x, y, z, pos1) :
    return sqrt((x-pos1[0])**2+(y-pos1[1])**2+(z-pos1[2])**2) < 100

           
  cMax=[-100000]*10
  cMin=[100000]*10
  
  n1_isIRPC = 0
  n2_isIRPC = 0
  
  Max=[-100000]*10
  Min=[100000]*10
  global HSCPmass
  HSCPmass = 0
  def timePlots (chain,HistL,Min,Max, Pileup = False):
    counter =0 

    for event in chain:

      n=len(event.simDigi1_isIRPC)
      n2=len(event.simDigi2_isIRPC)
      m=len(event.rpcHit_isIRPC)
      ngem = event.gemSegment_n

      pos11=(0,0,0)
      pos12=(0,0,0)
      pos21=(0,0,0)
      pos22=(0,0,0)

      r011,r012,r021,r022=[0,0,0,0]
      
      if m > 0 : counter = counter+1

      for i in range(n):
        if event.simDigi1_isIRPC[i]:
          if (abs(event.simDigi1_z[i])<1000 and abs(event.simDigi1_z[i])>900):
            pos11=(event.simDigi1_x[i],event.simDigi1_y[i],event.simDigi1_z[i])
            r011 = (event.simDigi1_tof[i]- event.simDigi1_t0[i])*SpeedOfLight
            # r011 can be fixed to 980 for preliminary noise study
#            r = sqrt(event.simDigi1_x[i]*event.simDigi1_x[i] + event.simDigi1_y[i]*event.simDigi1_y[i])
            z = event.simDigi1_z[i]
#            print z, " ", r011, 
            r = math.sqrt(r011*r011 - z*z)
            HistL[12].Fill(z,r)
          elif abs(event.simDigi1_z[i])>1000 :
            pos12=(event.simDigi1_x[i],event.simDigi1_y[i],event.simDigi1_z[i])
            r012 = (event.simDigi1_tof[i]- event.simDigi1_t0[i])*SpeedOfLight
            # r012 can be fixed to 1070 for preliminary noise study
#            z = event.simDigi1_z[i]
#            r = math.sqrt(r011*r011 - z*z)
           
#            HistL[12].Fill(z,r)
      for i in range(n2):
        if event.simDigi2_isIRPC[i]:
          if (abs(event.simDigi2_z[i])<1000 and abs(event.simDigi2_z[i])>900):
            pos21=(event.simDigi2_x[i],event.simDigi2_y[i],event.simDigi2_z[i])
            r021 = (event.simDigi2_tof[i]- event.simDigi2_t0[i])*SpeedOfLight
          elif abs(event.simDigi2_z[i])>1000 :
            pos22=(event.simDigi2_x[i],event.simDigi2_y[i],event.simDigi2_z[i])
            r022 = (event.simDigi2_tof[i]- event.simDigi2_t0[i])*SpeedOfLight

      t01,t02,tr1,tr2,dt,chi=[0,0,0,0,0,0]
      
      #geometry
      fakeRE31_plus = []
      fakeRE41_plus = []
      fakeRE31_minus = []
      fakeRE41_minus = []

      for i in range(n):
        if event.simDigi1_isIRPC[i]:
          if (abs(event.simDigi1_z[i])<1000 and abs(event.simDigi1_z[i])>900):
#            print r011
            HistL[10].Fill(r011/100.)
            t01=(event.simDigi1_t0[i])
          elif abs(event.simDigi1_z[i])>1000 :
#            print r012
            HistL[11].Fill(r012/100.)
            t02=(event.simDigi1_t0[i])

      RE31 = []
      RE41 = [] 

      for i in range(m):
        if event.rpcHit_isIRPC[i]:
          selection = False
          if Pileup:
               selection = (not isWithinTheCone(event.rpcHit_x[i],event.rpcHit_y[i],event.rpcHit_z[i],pos11) ) and (not isWithinTheCone(event.rpcHit_x[i],event.rpcHit_y[i],event.rpcHit_z[i],pos21) )
          else:
               selection = findClosestDigi(event.rpcHit_x[i], event.rpcHit_y[i],event.rpcHit_z[i],pos11,pos21)
          
          if abs(event.rpcHit_z[i])<1000 and abs(event.rpcHit_z[i])>900 and selection:
               tr1=event.rpcHit_time[i]-(math.sqrt(event.rpcHit_y[i]**2+event.rpcHit_x[i]**2+event.rpcHit_z[i]**2)-r011)/SpeedOfLight
               #using TVector3
               RE31.append( ROOT.TVector3( event.rpcHit_x[i], event.rpcHit_y[i], event.rpcHit_z[i] ) )
               #geometry
               if event.rpcHit_z[i] > 0: 
                 fakeRE31_plus.append(i)
               else : 
                 fakeRE31_minus.append(i)

          elif abs(event.rpcHit_z[i])>1000 and selection:
               #using TVector3
               RE41.append( ROOT.TVector3( event.rpcHit_x[i], event.rpcHit_y[i], event.rpcHit_z[i] ) )
               #geometry
               if event.rpcHit_z[i] > 0:
                 fakeRE41_plus.append(i)
               else : 
                 fakeRE41_minus.append(i) 

      ###### GEM association ######

      GE11 = []
      GE21 = []

      for i in range(ngem):
        if event.gemSegment_z[i] < 700:
          GE11.append( ROOT.TVector3( event.gemSegment_x[i], event.gemSegment_y[i], event.gemSegment_z[i] ) )
        else:
          GE21.append( ROOT.TVector3( event.gemSegment_x[i], event.gemSegment_y[i], event.gemSegment_z[i] ) )

      nRE31 = len(RE31)
      nRE41 = len(RE41)
      nGE11 = len(GE11)
      nGE21 = len(GE21)

      #print nRE41 , " = ", len(fakeRE41_plus), "+", len(fakeRE41_minus)
      #print "nRE31 = ", nRE31, " nRE41 = ", nRE41 , " nGE11 = " , nGE11, " nGE21 = " , nGE21

      dR_GE11_RE31 = 9
      dR_GE21_RE31 = 9

      for i in range(nRE31):
        for j in range(nGE11):
          tmp_dR = GE11[j].DeltaR(RE31[i])
          if tmp_dR < dR_GE11_RE31:
            dR_GE11_RE31 = tmp_dR
        for j in range(nGE21):
          tmp_dR = GE21[j].DeltaR(RE31[i])
          if tmp_dR < dR_GE21_RE31:
            dR_GE21_RE31 = tmp_dR

      #print "dR = " , dR_GE11_RE31, " " , dR_GE21_RE31

      HistL[12].Fill( dR_GE11_RE31 )
      HistL[13].Fill( dR_GE21_RE31 )

      ###############################

      l31_plus=len(fakeRE31_plus)
      l31_minus=len(fakeRE31_minus)

      l41_plus=len(fakeRE41_plus)
      l41_minus=len(fakeRE41_minus)

      phij = 0;
      phii = 0;
      etaj = 0;
      etai = 0;
      dphi = 0;
      deta = 0;

      noFilledChi = 1
      noFilledChi_select = 1

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
      
          dphi = phij-phii
          if dphi >= math.pi:
            dphi = dphi - 2*math.pi
          if dphi < (-1)*(math.pi):
            dphi = dphi + 2*math.pi

          dphi = dphi
          #dphi = RE41[0].DeltaPhi( RE31[0] )
          #dphi = abs(dphi)
          deta = etai-etaj

          if tr1 and tr2:
            dt=tr2-tr1
            chi=tr2**2+tr1**2
            if noFilledChi:
              HistL[5].Fill(chi)
              noFilledChi = 0
#deta cut -> if HSCP_noPU = 0.1, if HSCPLGW25 = 0.05
            if dt > -5 and abs(dphi) < 0.005 and abs(deta) < 0.05 and chi > 30:
              if noFilledChi_select:
                HistL[6].Fill(chi)
                noFilledChi_select = 0

      HistL[9].Fill(l31_plus*l41_plus)

      if t01:
        HistL[0].Fill(t01)
        if t01<Min[0]:Min[0]=t01
        elif t01>Max[0]:Max[0]=t01
      if t02:
        HistL[1].Fill(t02)
        if t02<Min[1]:Min[1]=t02
        elif t02>Max[1]:Max[1]=t02
      if tr1:
        HistL[2].Fill(tr1)
        if tr1<Min[2]:Min[2]=tr1
        elif tr1>Max[2]:Max[2]=tr1
      if tr2:
        HistL[3].Fill(tr2)
        if tr2<Min[3]:Min[3]=tr2
        elif tr2>Max[3]:Max[3]=tr2
'''
      if tr1 and tr2:
        dt=tr2-tr1
        HistL[4].Fill(dt)
        if dt<Min[4]:Min[4]=dt
        elif dt>Max[4]:Max[4]=dt
        chi=tr2**2+tr1**2
        HistL[5].Fill(chi)
        if dt > -5 and abs(dphi) < 0.005 and abs(deta) < 0.05:
          HistL[6].Fill(chi)
'''
      if dphi and deta and l31_plus > 0 and l41_plus > 0 and l31_plus == 1 and l41_plus == 1:
        HistL[7].Fill(dphi)
        HistL[8].Fill(deta)

    print counter   


  print "doing HSCP"
  timePlots(mychain,seriousHistL,Min,Max)
  print " "
  print "doing Muons"
  timePlots(mychain2,cseriousHistL,Min,Max, Pileup )
  
  print " All noise ", cseriousHistL[5].Integral(), " after cuts ", cseriousHistL[6].Integral()

  print (Min, Max)
  #print ("muons :",cMin, cMax)

  c1=ROOT.TCanvas()
  seriousHistL[0].SetLineColor(ROOT.kRed)
  gStyle.SetOptStat(0);
  seriousHistL[0].GetXaxis().SetTitle('Time delay')
  seriousHistL[0].GetYaxis().SetTitle('Number of Particles')
  seriousHistL[0].GetYaxis().SetTitleOffset(1.2)
  seriousHistL[0].SetTitle('First chamber, HSCP')
  seriousHistL[2].SetLineColor(ROOT.kBlue)
  seriousHistL[0].Draw("HIST")
  seriousHistL[2].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(seriousHistL[0],"Digital Hits","l")
  legend.AddEntry(seriousHistL[2],"RPC_Hits","l")
  legend.Draw() 
  c1.SaveAs(str_save+"timeHSCPFirst"+str_form)

  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[0].SetLineColor(ROOT.kRed)
  cseriousHistL[0].GetXaxis().SetTitle('Time delay')
  cseriousHistL[0].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[0].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[0].SetTitle('First chamber, muons')
  cseriousHistL[0].Draw("HIST")
  cseriousHistL[2].SetLineColor(ROOT.kBlue)
  cseriousHistL[2].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[0],"Digital Hits","l")
  legend.AddEntry(cseriousHistL[2],"RPC_Hits","l")
  legend.Draw()
#  Num_Muon1 = cseriousHist[]
  c1.SaveAs(str_save+"timeMuonsFirst"+str_form)

  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  seriousHistL[1].SetLineColor(ROOT.kRed)
  seriousHistL[1].GetXaxis().SetTitle('Time delay')
  seriousHistL[1].GetYaxis().SetTitle('Number of Particles')
  seriousHistL[1].GetYaxis().SetTitleOffset(1.2)
  seriousHistL[1].SetTitle('Second chamber, HSCP')
  seriousHistL[3].SetLineColor(ROOT.kBlue)
  seriousHistL[1].Draw("HIST")
  seriousHistL[3].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(seriousHistL[1],"Digital Hits","l")
  legend.AddEntry(seriousHistL[3],"RPC_Hits","l")
  legend.Draw() 
  c1.SaveAs(str_save+"timeHSCPSecond"+str_form)

  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[1].SetLineColor(ROOT.kRed)
  cseriousHistL[1].GetXaxis().SetTitle('Time delay')
  cseriousHistL[1].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[1].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[1].SetTitle('Second chamber, muons')
  cseriousHistL[1].Draw("HIST")
  cseriousHistL[3].SetLineColor(ROOT.kBlue)
  cseriousHistL[3].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[1],"Digital Hits","l")
  legend.AddEntry(cseriousHistL[3],"RPC_Hits","l")
  legend.Draw() 
  c1.SaveAs(str_save+"timeMuonsSecond"+str_form)

  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[4].SetLineColor(ROOT.kRed)
  cseriousHistL[4].SetLineStyle(2)
  cseriousHistL[4].GetXaxis().SetTitle('Time delay difference')
  cseriousHistL[4].GetXaxis().SetRangeUser(-20, 20)
  cseriousHistL[4].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[4].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[4].SetTitle('Difference between second and first chamber')
  seriousHistL[4].SetLineColor(ROOT.kBlue)
  seriousHistL[4].SetLineStyle(1)


  norm_mu = cseriousHistL[4].Integral()  # MG normalize to unity
  cseriousHistL[4].Scale(1./norm_mu)  # MG normalize to unity
  norm_hscp = seriousHistL[4].Integral()  # MG normalize to unity
  seriousHistL[4].Scale(1./norm_hscp)  # MG normalize to unity

  background_name = "Muon"
  if Pileup:
    background_name = "PileUp"
 
  cseriousHistL[4].Draw("HIST")
  seriousHistL[4].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[4],background_name,"l")
  legend.AddEntry(seriousHistL[4],"HSCP","l")
  legend.Draw()
  #line = ROOT.TLine(-5,0,-5,780)
  line = ROOT.TLine(-5,0,-5,0.56)
  line.Draw()

  #cseriousHistL[4].SetMaximum(0.3)

  c1.SaveAs(str_save+"timeDiff"+str_form)

  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[5].SetLineColor(ROOT.kRed)
  seriousHistL[5].SetLineColor(ROOT.kBlue)
  cseriousHistL[5].SetLineStyle(2)
  seriousHistL[5].SetLineStyle(1)
  cseriousHistL[5].GetXaxis().SetTitle('Chi2')
  cseriousHistL[5].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[5].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[5].SetTitle('Square sum of second and first chamber')

  print "Number of Muon before cut : ",cseriousHistL[5].GetEntries()

  norm_mu = cseriousHistL[5].Integral()  # MG normalize to unity
  cseriousHistL[5].Scale(1./norm_mu)  # MG normalize to unity
  norm_hscp = seriousHistL[5].Integral()  # MG normalize to unity
  seriousHistL[5].Scale(1./norm_hscp)  # MG normalize to unity



  cseriousHistL[5].Draw("HIST")
  seriousHistL[5].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[4],background_name,"l")
  legend.AddEntry(seriousHistL[4],"HSCP","l")
  legend.Draw() 
  c1.SaveAs(str_save+"timeChi"+str_form)


  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[6].SetLineColor(ROOT.kRed)
  cseriousHistL[6].GetXaxis().SetTitle('Chi2')
  cseriousHistL[6].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[6].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[6].SetTitle('Square sum of second and first chamber for positive propagation')
  cseriousHistL[6].SetLineStyle(2)
  seriousHistL[6].SetLineStyle(1)
  print "Number of Muon Aftere cut : ",cseriousHistL[6].GetEntries()
  cseriousHistL[6].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[6].Draw("HIST")

  seriousHistL[6].SetLineColor(ROOT.kBlue)
  seriousHistL[6].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[6].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[4],background_name,"l")
  legend.AddEntry(seriousHistL[4],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"timeChi_pos"+str_form)

  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[7].SetLineColor(ROOT.kRed)
  cseriousHistL[7].GetXaxis().SetTitle('dphi')
  cseriousHistL[7].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[7].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[7].SetTitle('Angle between second and first chamber for positive propagation')
  cseriousHistL[7].SetLineStyle(2)
  seriousHistL[7].SetLineStyle(1)
  cseriousHistL[7].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[7].Draw("HIST")
  cseriousHistL[7].SetMaximum(0.27)

  seriousHistL[7].SetLineColor(ROOT.kBlue)
  seriousHistL[7].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[7].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[7],background_name,"l")
  legend.AddEntry(seriousHistL[7],"HSCP","l")
  legend.Draw()
  line1 = ROOT.TLine(-0.005,0,-0.005,0.09)
  line1.Draw()
  line2 = ROOT.TLine(0.005,0,0.005,0.09)
  line2.Draw()
  c1.SaveAs(str_save+"dphi_pos_test2"+str_form)

  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[8].SetLineColor(ROOT.kRed)
  cseriousHistL[8].GetXaxis().SetTitle('deta')
  cseriousHistL[8].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[8].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[8].SetTitle('Angle between second and first chamber for positive propagation')
  cseriousHistL[8].SetLineStyle(2)
  seriousHistL[8].SetLineStyle(1)
  cseriousHistL[8].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[8].Draw("HIST")
  cseriousHistL[8].SetMaximum(0.27)

  seriousHistL[8].SetLineColor(ROOT.kBlue)
  seriousHistL[8].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[8].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[8],background_name,"l")
  legend.AddEntry(seriousHistL[8],"HSCP","l")
  legend.Draw()
  #line1 = ROOT.TLine(-0.05,0,-0.05,50)
  line1 = ROOT.TLine(-0.05,0,-0.05,0.12)
  line1.Draw()
  line2 = ROOT.TLine(0.05,0,0.05,0.12)
  #line2 = ROOT.TLine(0.05,0,0.05,50)
  line2.Draw()

  print seriousHistL[8].GetEntries()

  c1.SaveAs(str_save+"deta_pos"+str_form)

  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[9].SetLineColor(ROOT.kRed)
  cseriousHistL[9].GetXaxis().SetTitle('Number of Pairs')
  cseriousHistL[9].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[9].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[9].SetTitle('Npairs')
  cseriousHistL[9].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[9].Draw("HIST")
#  cseriousHistL[9].SetMaximum(1)


  seriousHistL[9].SetLineColor(ROOT.kBlue)
  seriousHistL[9].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[9].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[9],background_name,"l")
  legend.AddEntry(seriousHistL[9],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"npairs"+str_form)





  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[10].SetLineColor(ROOT.kRed)
  cseriousHistL[10].GetXaxis().SetTitle('Distance to the center of chamber(m)')
  cseriousHistL[10].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[10].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[10].SetTitle('Npairs')
#  cseriousHistL[10].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[10].Draw("HIST")
#  cseriousHistL[10].SetMaximum(1)


  seriousHistL[10].SetLineColor(ROOT.kBlue)
#  seriousHistL[10].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[10].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[10],background_name,"l")
  legend.AddEntry(seriousHistL[10],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"RE31"+str_form)


  c1=ROOT.TCanvas()
  cseriousHistL[12].SetLineColor(ROOT.kRed)
  cseriousHistL[12].SetLineStyle(2)
  cseriousHistL[12].GetXaxis().SetTitle('dR between GE11 and RE31')
  cseriousHistL[12].GetYaxis().SetTitle('Number of Particles')
#  cseriousHistL[12].SetTitle('Npairs')
  print "Number of Muon after GEM11 cut : ",cseriousHistL[12].GetEntries()
  cseriousHistL[12].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[12].Draw("HIST")
  cseriousHistL[12].SetMaximum(1)


  seriousHistL[12].SetLineColor(ROOT.kBlue)
  seriousHistL[12].SetLineStyle(1)
  seriousHistL[12].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[12].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[12],background_name,"l")
  legend.AddEntry(seriousHistL[12],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"dR_GE11_RE31"+str_form)

  c1=ROOT.TCanvas()
  cseriousHistL[13].SetLineColor(ROOT.kRed)
  cseriousHistL[13].SetLineStyle(2)
  cseriousHistL[13].GetXaxis().SetTitle('dR between GE21 and RE31')
  cseriousHistL[13].GetYaxis().SetTitle('Number of Particles')
#  cseriousHistL[13].SetTitle('Npairs')
  print "Number of Muon after GEM21 cut : ",cseriousHistL[13].GetEntries()
  cseriousHistL[13].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[13].Draw("HIST")
  cseriousHistL[13].SetMaximum(1)


  seriousHistL[13].SetLineColor(ROOT.kBlue)
  seriousHistL[13].SetLineStyle(1)
  seriousHistL[13].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[13].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[13],background_name,"l")
  legend.AddEntry(seriousHistL[13],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"dR_GE21_RE31"+str_form)

  c1=ROOT.TCanvas()
  cseriousHistL[14].Draw("COLZ")
#  cseriousHistL[10].SetMaximum(1)

  c1.SaveAs(str_save+"RECHAMBERS"+str_form)





  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0); 
  cseriousHistL[11].SetLineColor(ROOT.kRed)
  cseriousHistL[11].GetXaxis().SetTitle('Distance to the center of chamber(m)')
  cseriousHistL[11].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[11].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[11].SetTitle('Npairs')
#  cseriousHistL[11].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[11].Draw("HIST")
#  cseriousHistL[11].SetMaximum(1)


  seriousHistL[11].SetLineColor(ROOT.kBlue)
#  seriousHistL[11].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[11].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[11],background_name,"l")
  legend.AddEntry(seriousHistL[11],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"RE41"+str_form)



  c1=ROOT.TCanvas()
  propHist=ROOT.TGraph(numBin[5])
  n1 = seriousHistL[5].Integral()
  n2 = cseriousHistL[5].Integral()
 
  count1=n1
  count2=n2
  print("nentries hist5 :",n1,n2)
  for q in range(numBin[5]):
      count1-=seriousHistL[5].GetBinContent(q)
      count2-=cseriousHistL[5].GetBinContent(q)
      propHist.SetPoint(q,count1/n1,count2/n2)
      print(count1,count2,n1,n2)

# MG beggin: chi2 avec condition positive

  propHist_pos=ROOT.TGraph(numBin[6])
  n1_pos = seriousHistL[6].Integral()
  n2_pos = cseriousHistL[6].Integral()

  count1=n1_pos
  count2=n2_pos

  print("nentries hist6 :",n1_pos,n2_pos)
  for q in range(numBin[6]):
      count1-=seriousHistL[6].GetBinContent(q)
      count2-=cseriousHistL[6].GetBinContent(q)
      propHist_pos.SetPoint(q,count1/n1,count2/n2)
##  print(count1,count2,n1,n2)

# MG end

  propHist.SetTitle( 'Muons vs '+str(HSCPmass)+' HSCP kept')
  propHist.Draw()
  gStyle.SetOptStat(0);
  propHist.GetXaxis().SetTitle('HSCP')
  propHist.GetXaxis().SetRangeUser(0.,1.)
  propHist.GetYaxis().SetTitle(background_name)
  propHist.GetYaxis().SetTitleOffset(1.2)
  propHist.SetMarkerColor(ROOT.kBlue)
  propHist.SetMarkerStyle(21)
  propHist.SetLineStyle(2)
  propHist.Draw("PE")

  propHist.SetMinimum(1e-3)



  propHist_pos.SetMarkerColor(ROOT.kRed)
  propHist_pos.SetMarkerStyle(20)
  propHist_pos.SetLineStyle(1)
  propHist_pos.Draw("SAMEPE")


  legend = ROOT.TLegend(0.13,0.70,0.55,0.87)
  legend.AddEntry(propHist,"chi2 discri","p")
  legend.AddEntry(propHist_pos,"chi2 + (dT > -5 ns & |dphi| < 0.005 & |deta| < 0.05)","p")
  legend.Draw() 

  c1.SetLogy()
  c1.SaveAs(str_save+"graal"+str_form)


  
Pileup = False
readLorentzVector(Pileup)

