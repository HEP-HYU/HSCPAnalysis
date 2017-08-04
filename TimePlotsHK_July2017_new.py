#!/usr/bin/env python

import copy
import os, sys, ROOT, math
from math import exp, atan, sqrt, log, atan2, tan
from array import array
from ROOT import gStyle
#from rootpy.interactive import wait

Pileup = False
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
  str_save="plot_HSCPLGW25_noise_BgPU0/"
  #str_save="plot_HSCPLGW25_noise_BgSNBMAPU200/"
  #str_save = "plot_NOnoise_noPU/"
  #str_save="plots/new/Go"
  #str_save="plots/HSCP_old_NoPU/"
  #str_save="plots/910pre3/"
  #str_form=".root"
  str_form=".pdf"
  
  c1=ROOT.TCanvas()
  seriousHistL=[]
  cseriousHistL=[]

  seriousHistL.append(ROOT.TH1F("seriousHist0"  ,"placeholder", 50, -10, 100))  
  seriousHistL.append(ROOT.TH1F("seriousHist1"  ,"placeholder", 50, -10, 100))  
  seriousHistL.append(ROOT.TH1F("seriousHist2"  ,"placeholder", 50, -10, 100))  
  seriousHistL.append(ROOT.TH1F("seriousHist3"  ,"placeholder", 50, -10, 100))  
  seriousHistL.append(ROOT.TH1F("seriousHist4"  ,"placeholder", 50, -35, 85))  
  seriousHistL.append(ROOT.TH1F("seriousHist5"  ,"placeholder", 50, -10, 200))  
  seriousHistL.append(ROOT.TH1F("seriousHist6"  ,"placeholder", 50, -10, 200))  
  seriousHistL.append(ROOT.TH1F("seriousHist7"  ,"placeholder", 100, -0.011, 0.011))  
  seriousHistL.append(ROOT.TH1F("seriousHist8"  ,"placeholder", 100, -0.2, 0.2))  
  seriousHistL.append(ROOT.TH1F("seriousHist9"  ,"placeholder", 50, 0, 10))  
  seriousHistL.append(ROOT.TH1F("seriousHist10"  ,"placeholder", 50, 8, 13))  
  seriousHistL.append(ROOT.TH1F("seriousHist11"  ,"placeholder", 50, 8, 13))  
  seriousHistL.append(ROOT.TH1F("seriousHist12"  ,"placeholder", 100, 0, 3))  
  seriousHistL.append(ROOT.TH1F("seriousHist13"  ,"placeholder", 100, 0, 3))  
  seriousHistL.append(ROOT.TH2F("seriousHist14","placeholder",200,900,1100, 200, 200, 400))
  seriousHistL.append(ROOT.TH2F("seriousHist15","placeholder",200, 200, 400, 200, 200, 400))
  seriousHistL.append(ROOT.TH1F("seriousHist16","dR between rechit and simdigi1",100, 0, 1))
  seriousHistL.append(ROOT.TH1F("seriousHist17","dR between rechit and simdigi2",100, 0, 1))

  cseriousHistL.append(ROOT.TH1F("cseriousHist0"  ,"placeholder", 50, -10, 100))
  cseriousHistL.append(ROOT.TH1F("cseriousHist1"  ,"placeholder", 50, -10, 100))
  cseriousHistL.append(ROOT.TH1F("cseriousHist2"  ,"placeholder", 50, -10, 100))
  cseriousHistL.append(ROOT.TH1F("cseriousHist3"  ,"placeholder", 50, -10, 100))
  cseriousHistL.append(ROOT.TH1F("cseriousHist4"  ,"placeholder", 50, -35, 85))
  cseriousHistL.append(ROOT.TH1F("cseriousHist5"  ,"placeholder", 50, -10, 200))
  cseriousHistL.append(ROOT.TH1F("cseriousHist6"  ,"placeholder", 50, -10, 200))
  cseriousHistL.append(ROOT.TH1F("cseriousHist7"  ,"placeholder", 100, -0.011, 0.011))
  cseriousHistL.append(ROOT.TH1F("cseriousHist8"  ,"placeholder", 100, -0.2, 0.2))
  cseriousHistL.append(ROOT.TH1F("cseriousHist9"  ,"placeholder", 50, 0, 10))
  cseriousHistL.append(ROOT.TH1F("cseriousHist10"  ,"placeholder", 50, 8, 13))
  cseriousHistL.append(ROOT.TH1F("cseriousHist11"  ,"placeholder", 50, 8, 13))
  cseriousHistL.append(ROOT.TH1F("cseriousHist12"  ,"placeholder", 100, 0, 3))
  cseriousHistL.append(ROOT.TH1F("cseriousHist13"  ,"placeholder", 100, 0, 3))
  cseriousHistL.append(ROOT.TH2F("cseriousHist14","placeholder",200,900,1100, 200, 200, 400))
  cseriousHistL.append(ROOT.TH2F("cseriousHist15","placeholder",200, 200, 400, 200, 200, 400))
  cseriousHistL.append(ROOT.TH1F("cseriousHist16","dR between rechit and simdigi1",100, 0, 1))
  cseriousHistL.append(ROOT.TH1F("cseriousHist17","dR between rechit and simdigi2",100, 0, 1))
 
  def findClosestDigi (x, y, z, pos1, pos2) :
    return (x-pos1[0])**2+(y-pos1[1])**2+(z-pos1[2])**2<(x-pos2[0])**2+(y-pos2[1])**2+(z-pos2[2])**2

  def isWithinTheCone (x, y, z, pos1) :
    return sqrt((x-pos1[0])**2+(y-pos1[1])**2+(z-pos1[2])**2) < 100

  def DrawCanvas( h1, cor1, l1, h2, cor2, l2, xtitle, ytitle, name, norm = True, PlotMin, PlotMax, LineX1, LineY1, LineX2, LineY2) :

###h1 -> HSCP Histogram, h2 -> Muon(or PU) Histogram
     c1=ROOT.TCanvas()
     h1.SetLineColor(cor1)
     gStyle.SetOptStat(0);
     h2.GetXaxis().SetTitle(xtitle)
     h2.GetYaxis().SetTitle(ytitle)
     h2.GetYaxis().SetTitleOffset(1.2)
     h2.SetTitle('')
     h2.SetLineColor(cor2)
     h1.SetLineStyle(1)
     h2.SetLineStyle(2)
     h1.SetLineWidth(2)
     h2.SetLineWidth(2)
     if norm:
       if not h1.Integral() == 0 :
         h1.Scale(1.0/h1.Integral())
       if not h2.Integral() == 0 :
         h2.Scale(1.0/h2.Integral())
     h2.Draw("HIST") #Muon(orPU)
     if not PlotMin  -1.0 :
       h2.SetMinimum(PlotMin)
     if not PlotMax == -1.0 :
       h2.SetMaximum(PlotMax)
     legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
     if Pileup == False :
       h1.Draw("SAME HIST") #HSCP
       legend.AddEntry(h1,l1,"l")
     legend.AddEntry(h2,l2,"l")
     legend.SetLineColor(0);
     legend.SetFillColor(0);
     legend.Draw()
     line1 = ROOT.TLine(LineX1, LineY1, LineX2, LineY2)
     line2 = ROOT.TLine(-LineX1, LineY1, -LineX2, LineY2)
     line1.Draw()
     if (h1 == seriousHistL[7] or h1 == seriousHistL[8]) :
       line2.Draw()
#     c1.SaveAs(name)
     c1.SaveAs(str_save+name+str_form)

           
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
    counter_pre = 0
    counter_step1 = 0
    counter_step2 = 0
    counter_step3 = 0

    for event in chain:

      n=len(event.simDigi1_isIRPC)
      n2=len(event.simDigi2_isIRPC)
      m=len(event.rpcHit_isIRPC)
      ngem = event.gemSegment_n

      simDigi1_RE31 = ROOT.TLorentzVector(0,0,0,0)
      simDigi1_RE41 = ROOT.TLorentzVector(0,0,0,0)
      simDigi2_RE31 = ROOT.TLorentzVector(0,0,0,0)
      simDigi2_RE41 = ROOT.TLorentzVector(0,0,0,0)
     
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
            simDigi1_RE31 = ROOT.TLorentzVector( event.simDigi1_x[i],event.simDigi1_y[i],event.simDigi1_z[i], event.simDigi1_tof[i] ) 
            r011 = (event.simDigi1_tof[i]- event.simDigi1_t0[i])*SpeedOfLight

            HistL[10].Fill(r011/100.)
            t01=(event.simDigi1_t0[i])

            # r011 can be fixed to 980 for preliminary noise study
#            r = sqrt(event.simDigi1_x[i]*event.simDigi1_x[i] + event.simDigi1_y[i]*event.simDigi1_y[i])
            z = event.simDigi1_z[i]
#            print z, " ", r011, 
            r = math.sqrt(r011*r011 - z*z)
            HistL[12].Fill(z,r)
          elif abs(event.simDigi1_z[i])>1000 :
            pos12=(event.simDigi1_x[i],event.simDigi1_y[i],event.simDigi1_z[i])
            simDigi1_RE41 = ROOT.TLorentzVector( event.simDigi1_x[i],event.simDigi1_y[i],event.simDigi1_z[i], event.simDigi1_tof[i] )
            r012 = (event.simDigi1_tof[i]- event.simDigi1_t0[i])*SpeedOfLight

            HistL[11].Fill(r012/100.)
            t02=(event.simDigi1_t0[i])

            # r012 can be fixed to 1070 for preliminary noise study
#            z = event.simDigi1_z[i]
#            r = math.sqrt(r011*r011 - z*z)
           
#            HistL[12].Fill(z,r)
      for i in range(n2):
        if event.simDigi2_isIRPC[i]:
          if (abs(event.simDigi2_z[i])<1000 and abs(event.simDigi2_z[i])>900):
            pos21=(event.simDigi2_x[i],event.simDigi2_y[i],event.simDigi2_z[i])
            simDigi2_RE31 = ROOT.TLorentzVector( event.simDigi2_x[i],event.simDigi2_y[i],event.simDigi2_z[i], event.simDigi2_tof[i] )
            r021 = (event.simDigi2_tof[i]- event.simDigi2_t0[i])*SpeedOfLight
          elif abs(event.simDigi2_z[i])>1000 :
            pos22=(event.simDigi2_x[i],event.simDigi2_y[i],event.simDigi2_z[i])
            simDigi2_RE41 = ROOT.TLorentzVector( event.simDigi2_x[i],event.simDigi2_y[i],event.simDigi2_z[i], event.simDigi2_tof[i] )
            r022 = (event.simDigi2_tof[i]- event.simDigi2_t0[i])*SpeedOfLight

      t01,t02,tr1,tr2,dt,chi=[0,0,0,0,0,0]
     
      RE31 = []
      RE41 = [] 

      #### RPC hits ####
      for i in range(m):
        if event.rpcHit_isIRPC[i]:
          tmpLorentzVector = ROOT.TLorentzVector( event.rpcHit_x[i], event.rpcHit_y[i], event.rpcHit_z[i], event.rpcHit_time[i])

          dR_digi1_RE31 = tmpLorentzVector.DeltaR( simDigi1_RE31 )
          dR_digi2_RE31 = tmpLorentzVector.DeltaR( simDigi2_RE31 )
          dR_digi1_RE41 = tmpLorentzVector.DeltaR( simDigi1_RE41 )
          dR_digi2_RE41 = tmpLorentzVector.DeltaR( simDigi2_RE41 )

          selection = False
          if Pileup:
               selection = (not isWithinTheCone(event.rpcHit_x[i],event.rpcHit_y[i],event.rpcHit_z[i],pos11) ) and (not isWithinTheCone(event.rpcHit_x[i],event.rpcHit_y[i],event.rpcHit_z[i],pos21) )
          else:
               if abs(event.rpcHit_z[i])<1000 and abs(event.rpcHit_z[i])>900:
                 HistL[16].Fill(dR_digi1_RE31)
                 HistL[17].Fill(dR_digi2_RE31)
               elif abs(event.rpcHit_z[i])>1000:
                 HistL[16].Fill(dR_digi1_RE41)  
                 HistL[17].Fill(dR_digi2_RE41)  
               selection = dR_digi1_RE31 < 0.1 or dR_digi2_RE31 < 0.1 or dR_digi1_RE41 < 0.1 or dR_digi2_RE41 < 0.1
               #selection = findClosestDigi(event.rpcHit_x[i], event.rpcHit_y[i],event.rpcHit_z[i],pos11,pos21) or findClosestDigi(event.rpcHit_x[i], event.rpcHit_y[i],event.rpcHit_z[i],pos12,pos22)
          
          if abs(event.rpcHit_z[i])<1000 and abs(event.rpcHit_z[i])>900 and selection:
               RE31.append( ROOT.TLorentzVector( event.rpcHit_x[i], event.rpcHit_y[i], event.rpcHit_z[i], event.rpcHit_time[i]) )
          elif abs(event.rpcHit_z[i])>1000 and selection:
               RE41.append( ROOT.TLorentzVector( event.rpcHit_x[i], event.rpcHit_y[i], event.rpcHit_z[i], event.rpcHit_time[i] ) )

      ###### GEM association ######

      GE11 = []
      GE21 = []

      for i in range(ngem):
        if event.gemSegment_z[i] < 700:
          GE11.append( ROOT.TLorentzVector( event.gemSegment_x[i], event.gemSegment_y[i], event.gemSegment_z[i], event.gemSegment_time[i] ) )
        else:
          GE21.append( ROOT.TLorentzVector( event.gemSegment_x[i], event.gemSegment_y[i], event.gemSegment_z[i], event.gemSegment_time[i] ) )

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

      FoundPair = False

      for i in range(nRE31):
        for j in range(nRE41):
          tr1 = RE31[i].T() - (RE31[i].P()-r012)/SpeedOfLight
          tr2 = RE41[j].T() - (RE41[j].P()-r012)/SpeedOfLight
          dphi = RE41[0].DeltaPhi( RE31[0] )
          deta = RE31[0].Eta() - RE41[0].Eta() 
 
          dt=tr2-tr1
          HistL[4].Fill(dt)

          chi=tr2**2+tr1**2

          HistL[5].Fill(chi)
          HistL[7].Fill(dphi)
          HistL[8].Fill(deta)

          counter_pre = counter_pre + 1 
          step1 = abs(dphi) < 0.005 and abs(deta) < 0.05 and chi > 30 #deta cut -> if HSCP_noPU = 0.1, if HSCPLGW25 = 0.05
          step2 = dR_GE11_RE31 < 0.2 and dR_GE21_RE31 and step1
          if dR_GE11_RE31 > 8:
            step2 = dR_GE21_RE31 and step1
          step3 = dt > -5 and step2

          if step1:
            counter_step1 = counter_step1 + 1

          if step2:
            counter_step2 = counter_step2 + 1

          if step3: 
            counter_step3 = counter_step3 + 1
            HistL[6].Fill(chi)

          FoundPair = True
          break
        if FoundPair:
          break

      HistL[9].Fill(nRE31*nRE41)

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

#      if tr1 and tr2:
#        dt=tr2-tr1
#        HistL[4].Fill(dt)
#        if dt<Min[4]:Min[4]=dt
#        elif dt>Max[4]:Max[4]=dt
#        chi=tr2**2+tr1**2
#        HistL[5].Fill(chi)
#        if dt > -5 and abs(dphi) < 0.005 and abs(deta) < 0.05:
#          HistL[6].Fill(chi)

    print "total = ", counter   
    print "foundpair = ", counter_pre
    print "step1 = ", counter_step1
    print "step2 = ", counter_step2
    print "step3 = ", counter_step3
   


  print "doing HSCP"
  timePlots(mychain,seriousHistL,Min,Max)
  print " "
  print "doing Muons"
  timePlots(mychain2,cseriousHistL,Min,Max, Pileup )
  
#  print " All noise ", cseriousHistL[5].Integral(), " after cuts ", cseriousHistL[6].Integral()

  print (Min, Max)
  #print ("muons :",cMin, cMax)
'''
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
'''
  DrawCanvas(seriousHistL[2], ROOT.kBlue, "RPC_Hits", seriousHistL[0], ROOT.kRed, "Digital Hits", 'Time delay', 'Number of Particles', "timeHSCPFirst", False, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)
  
'''
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
'''

  DrawCanvas(cseriousHistL[2], ROOT.kBlue, "RPC_Hits", cseriousHistL[0], ROOT.kRed, "Digital Hits", 'Time delay', 'Number of Particles', "timeMuonsFirst", False, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)

'''
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
'''

  DrawCanvas(seriousHistL[3], ROOT.kBlue, "RPC_Hits", seriousHistL[1], ROOT.kRed, "Digital Hits", 'Time delay', 'Number of Particles', "timeHSCPSecond", False, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)

'''
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
'''

  DrawCanvas(cseriousHistL[3], ROOT.kBlue, "RPC_Hits", cseriousHistL[1], ROOT.kRed, "Digital Hits", 'Time delay', 'Number of Particles', "timeMuonsSecond", False, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)

  background_name = "Muon"
  if Pileup:
    background_name = "PileUp"

'''
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
'''  
  print "timeDiff >> Number of ", background_name, " : ", cseriousHistL[4].Integral()
  print "timeDiff >> Number of HSCP : ", seriousHistL[4].Integral()
'''
  norm_mu = cseriousHistL[4].Integral()  # MG normalize to unity
  cseriousHistL[4].Scale(1./norm_mu)  # MG normalize to unity
  norm_hscp = seriousHistL[4].Integral()  # MG normalize to unity
  seriousHistL[4].Scale(1./norm_hscp)  # MG normalize to unity

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
'''
  if Pileup == False :
    timeDiff_Y2 = 0.56
  else :
    timeDiff_Y2 = 0.14
  DrawCanvas(seriousHistL[4], ROOT.kBlue, "HSCP", cseriousHistL[4], ROOT.kRed, background_name, 'Time delay difference', 'Number of Particles', "timeDiff", True, 0.0, 0.15, -5.0, 0.0, -5.0, timeDiff_Y2)

'''
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
'''
  print "chi2 >> Number of ", background_name," before cut : ",cseriousHistL[5].Integral()
  print "chi2 >> Number of HSCP before cut : ", seriousHistL[5].Integral()
'''
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
'''
  DrawCanvas(seriousHistL[5], ROOT.kBlue, "HSCP", cseriousHistL[5], ROOT.kRed, background_name, r'$\Chi^2$', 'Number of Particles', "timeChi", True, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)
'''
  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[6].SetLineColor(ROOT.kRed)
  cseriousHistL[6].GetXaxis().SetTitle('Chi2')
  cseriousHistL[6].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[6].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[6].SetTitle('Square sum of second and first chamber for positive propagation')
  cseriousHistL[6].SetLineStyle(2)
  seriousHistL[6].SetLineStyle(1)
'''  
  print "chi2 >> Number of", background_name, " Aftere cut(chi2>30) : ",cseriousHistL[6].Integral()
  print "chi2 >> Number of HSCP  Aftere cut(chi2>30) : ",seriousHistL[6].Integral()
'''
  norm_mu = cseriousHistL[6].Integral()
  cseriousHistL[6].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[6].Draw("HIST")
  cseriousHistL[6].SetMaximum(0.06)

  seriousHistL[6].SetLineColor(ROOT.kBlue)
  norm_hscp = seriousHistL[6].Integral()
  seriousHistL[6].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[6].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[4],background_name,"l")
  legend.AddEntry(seriousHistL[4],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"timeChi_pos"+str_form)
'''
  DrawCanvas(seriousHistL[6], ROOT.kBlue, "HSCP", cseriousHistL[6], ROOT.kRed, background_name, r'$\Chi^2$', 'Number of Particles', "timeChi_pos", True, -1.0, 0.06, 0.0, 0.0, 0.0, 0.0)
'''
  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[7].SetLineColor(ROOT.kRed)
  cseriousHistL[7].GetXaxis().SetTitle('dphi')
  cseriousHistL[7].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[7].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[7].SetTitle('Angle between second and first chamber for positive propagation')
  cseriousHistL[7].SetLineStyle(2)
  seriousHistL[7].SetLineStyle(1)
'''
  print "dphi >> Number of ",background_name," : ", cseriousHistL[7].Integral(28, 72)
  print "dphi >> Number of HSCP : ", seriousHistL[7].Integral(28, 72)
'''
  norm_mu = cseriousHistL[7].Integral()
  cseriousHistL[7].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[7].Draw("HIST")
  cseriousHistL[7].SetMaximum(0.27)

  seriousHistL[7].SetLineColor(ROOT.kBlue)
  norm_hscp = seriousHistL[7].Integral()
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
'''
  DrawCanvas(seriousHistL[7], ROOT.kBlue, "HSCP", cseriousHistL[7], ROOT.kRed, background_name, r'$\delta\phi$', 'Number of Particles', "dphi_pos", True, -1.0, 0.17, 0.005, 0.0, 0.005, 0.09)
'''
  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[8].SetLineColor(ROOT.kRed)
  cseriousHistL[8].GetXaxis().SetTitle('deta')
  cseriousHistL[8].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[8].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[8].SetTitle('Angle between second and first chamber for positive propagation')
  cseriousHistL[8].SetLineStyle(2)
  seriousHistL[8].SetLineStyle(1)
'''
  print "deta >> Number of ",background_name," : ", cseriousHistL[8].Integral(38, 62)
  print "deta >> Number of HSCP : ", seriousHistL[8].Integral(38, 62)
'''
  norm_mu = cseriousHistL[8].Integral()
  norm_hscp = seriousHistL[8].Integral()
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

  c1.SaveAs(str_save+"deta_pos"+str_form)
'''
  if Pileup == False :
    deta_Y2 = 0.17
  else :
    deta_Y2 = 0.03
  DrawCanvas(seriousHistL[8], ROOT.kBlue, "HSCP", cseriousHistL[8], ROOT.kRed, background_name, r'$\delta\eta$', 'Number of Particles', "deta_pos", True, -1.0, deta_Y2, 0.05, 0.0, 0.05, deta_Y2)
'''
  c1=ROOT.TCanvas()
  gStyle.SetOptStat(0);
  cseriousHistL[9].SetLineColor(ROOT.kRed)
  cseriousHistL[9].GetXaxis().SetTitle('Number of Pairs')
  cseriousHistL[9].GetYaxis().SetTitle('Number of Particles')
  cseriousHistL[9].GetYaxis().SetTitleOffset(1.2)
  cseriousHistL[9].SetTitle('Npairs')
'''
  print "Npairs >> Number of ",background_name," : ", cseriousHistL[9].Integral()
  print "Npairs >> Number of HSCP : ", seriousHistL[9].Integral()
'''
  norm_mu = cseriousHistL[9].Integral()
  norm_hscp = seriousHistL[9].Integral()
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
'''
  DrawCanvas(seriousHistL[9], ROOT.kBlue, "HSCP", cseriousHistL[9], ROOT.kRed, background_name, 'Number of pairs', 'Number of Particles', "npairs", True, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)
'''
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
'''
  DrawCanvas(seriousHistL[10], ROOT.kBlue, "HSCP", cseriousHistL[10], ROOT.kRed, background_name, 'Distance to the center of chamber(m)', 'Number of Particles', "RE31", False, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)
'''
  c1=ROOT.TCanvas()
  cseriousHistL[12].SetLineColor(ROOT.kRed)
  cseriousHistL[12].SetLineStyle(2)
  cseriousHistL[12].GetXaxis().SetTitle('dR between GE11 and RE31')
  cseriousHistL[12].GetYaxis().SetTitle('Number of Particles')
#  cseriousHistL[12].SetTitle('Npairs')
'''
  print "GEM11 >> Number of ",background_name,"after GEM11 cut : ", cseriousHistL[12].Integral()
  print "GEM11 >> Number of HSCP after GEM11 cut : ", seriousHistL[12].Integral()
'''
  norm_mu = cseriousHistL[12].Integral()
  norm_hscp = seriousHistL[12].Integral()
  cseriousHistL[12].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[12].Draw("HIST")
  cseriousHistL[12].SetMaximum(1)


  seriousHistL[12].SetLineColor(ROOT.kBlue)
  seriousHistL[12].SetLineStyle(1)
  if not norm_hscp == 0:
    seriousHistL[12].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[12].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[12],background_name,"l")
  legend.AddEntry(seriousHistL[12],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"dR_GE11_RE31"+str_form)
'''
  if Pileup == False :
    GE11_Y2 = 0.7
  else :
    GE11_Y2 = 0.28
  DrawCanvas(seriousHistL[12], ROOT.kBlue, "HSCP", cseriousHistL[12], ROOT.kRed, background_name, 'dR between GE11 and RE31', 'Number of Particles', "dR_GE11_RE31", True, -1.0, GE11_Y2, 0.0, 0.0, 0.0, 0.0)
'''
  c1=ROOT.TCanvas()
  cseriousHistL[13].SetLineColor(ROOT.kRed)
  cseriousHistL[13].SetLineStyle(2)
  cseriousHistL[13].GetXaxis().SetTitle('dR between GE21 and RE31')
  cseriousHistL[13].GetYaxis().SetTitle('Number of Particles')
#  cseriousHistL[13].SetTitle('Npairs')
'''
  print "GEM21 >> Number of ",background_name,"after GEM21 cut : ", cseriousHistL[13].Integral()
  print "GEM21 >> Number of HSCP after GEM21 cut : ", seriousHistL[13].Integral()
'''
  norm_mu = cseriousHistL[13].Integral()
  norm_hscp = seriousHistL[13].Integral()
  cseriousHistL[13].Scale(1./norm_mu)  # MG normalize to unity
  cseriousHistL[13].Draw("HIST")
  cseriousHistL[13].SetMaximum(1)


  seriousHistL[13].SetLineColor(ROOT.kBlue)
  seriousHistL[13].SetLineStyle(1)
  if not norm_hscp == 0:
    seriousHistL[13].Scale(1./norm_hscp)  # MG normalize to unity
  seriousHistL[13].Draw("SAME HIST")
  legend = ROOT.TLegend(0.67,0.70,0.87,0.87)
  legend.AddEntry(cseriousHistL[13],background_name,"l")
  legend.AddEntry(seriousHistL[13],"HSCP","l")
  legend.Draw()
  c1.SaveAs(str_save+"dR_GE21_RE31"+str_form)
'''
  if Pileup == False :
    GE21_Y2 = 0.6
  else :
    GE21_Y2 = 0.2
  DrawCanvas(seriousHistL[13], ROOT.kBlue, "HSCP", cseriousHistL[13], ROOT.kRed, background_name, 'dR between GE21 and RE31', 'Number of Particles', "dR_GE21_RE31", True, -1.0, GE11_Y2, 0.0, 0.0, 0.0, 0.0)


  c1=ROOT.TCanvas()
  cseriousHistL[14].Draw("COLZ")
#  cseriousHistL[10].SetMaximum(1)

  c1.SaveAs(str_save+"RECHAMBERS"+str_form)


'''
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
'''
  DrawCanvas(seriousHistL[11], ROOT.kBlue, "HSCP", cseriousHistL[11], ROOT.kRed, background_name, 'Distance to the center of chamber(m)', 'Number of Particles', "RE41", False, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)

  DrawCanvas( seriousHistL[16], ROOT.kBlue, "HSCP", cseriousHistL[16], ROOT.kRed, "Muons", "dR between rechit and simdigi1", "Entries", "dRsimdigi1", True, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)
  DrawCanvas( seriousHistL[17], ROOT.kBlue, "HSCP", cseriousHistL[17], ROOT.kRed, "Muons", "dR between rechit and simdigi2", "Entries", "dRsimdigi2", True, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0)

  c1=ROOT.TCanvas()
  propHist=ROOT.TGraph( seriousHistL[5].GetNbinsX() )
  n1 = seriousHistL[5].Integral()
  n2 = cseriousHistL[5].Integral()
 
  count1=n1
  count2=n2
  print("nentries hist5 :",n1,n2)
  for q in range( seriousHistL[5].GetNbinsX() ):
      count1-=seriousHistL[5].GetBinContent(q)
      count2-=cseriousHistL[5].GetBinContent(q)
      propHist.SetPoint(q,count1/n1,count2/n2)
      print(count1,count2,n1,n2)

# MG beggin: chi2 avec condition positive

  propHist_pos=ROOT.TGraph( seriousHistL[6].GetNbinsX()   )
  n1_pos = seriousHistL[6].Integral()
  n2_pos = cseriousHistL[6].Integral()

  count1=n1_pos
  count2=n2_pos

  print("nentries hist6 :",n1_pos,n2_pos)
  for q in range( seriousHistL[6].GetNbinsX()   ):
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


  legend = ROOT.TLegend(0.13,0.70,0.60,0.87)
  legend.AddEntry(propHist,"chi2 discri","p")
  legend.AddEntry(propHist_pos,"chi2 + (dT > -5 ns & |dphi| < 0.005 & |deta| < 0.05) & chi2 > 30","p")
  legend.Draw() 

  c1.SetLogy()
  c1.SaveAs(str_save+"graal"+str_form)


  
readLorentzVector(Pileup)

