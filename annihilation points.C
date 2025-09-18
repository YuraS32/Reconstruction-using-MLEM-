#define medicalProject_cxx
#include "medicalProject.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <math.h>
#include<TGraph2D.h>
#include<TH3D.h>

void medicalProject::Loop()
{
   if (fChain == 0) return;

   TFile fNew("output.root","RECREATE");


   TH1D *ann_x = new TH1D("ann_x","X of annihilation",90,-45,45); 
   TH1D *ann_y = new TH1D("ann_y","Y of annihilation",90,-45,45); 
   TH1D *ann_z = new TH1D("ann_z","Z of annihilation",90,-45,45); 

   TGraph2D *pl_3d = new TGraph2D();
   
   const Float_t c = 29.979;  // cm/ns;
    TH3D *h3d = new TH3D("h3d", "XYZ density",
                        60, -40, 40,
                        60, -40, 40,
                        60, -25, 25);

   TH2D *h_xy = new TH2D("h_xy", "XY density;X;Y",
                         100, -15, 15, 100, -15, 15);
   TH2D *h_yz = new TH2D("h_yz", "YZ density;Y;Z",
                         150, -40, 40, 150, -25, 25);
   TH2D *h_xz = new TH2D("h_xz", "XZ density;X;Z",
                         150, -40, 40, 150, -25, 25); 
    
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   int i=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
   
	   Long64_t ientry = LoadTree(jentry);
      
	   if (ientry < 0) break;
      
	   nb = fChain->GetEntry(jentry);   nbytes += nb;

	  // if(X_hit1==0 || X_hit2==0) continue;
       
      Float_t dx = X_hit2 - X_hit1;
      Float_t dy = Y_hit2 - Y_hit1;
      Float_t dz = Z_hit2 - Z_hit1;

      Float_t L = sqrt(dx*dx + dy*dy + dz*dz);
      if (L==0) continue;

      Float_t dt = t_hit2 - t_hit1;
      Float_t dL = c * dt / 2;

      Float_t ux = dx / L;
      Float_t uy = dy / L;
      Float_t uz = dz / L;

      Float_t X_ann = (X_hit1+X_hit2)/2 - ux * dL;
      Float_t Y_ann = (Y_hit1+Y_hit2)/2 - uy * dL;
      Float_t Z_ann = (Z_hit1+Z_hit2)/2 - uz * dL;
      i++;
      if ((X_ann<=-38.407) || (X_ann>=38.407) || (Y_ann<=-38.407) || (Y_ann>=38.407) || (Z_ann<=-23) || (Z_ann>=23))
      {
         continue;
      }
      ann_x->Fill(X_ann);
      ann_y->Fill(Y_ann);
      ann_z->Fill(Z_ann);
      if (i%1==0)
      {
          h3d->Fill(X_ann, Y_ann, Z_ann);
          h_xy->Fill(X_ann, Y_ann);
          h_yz->Fill(Y_ann, Z_ann);
          h_xz->Fill(X_ann, Z_ann);
      }

   }
      
      ann_x->Write();
      ann_y->Write();
      ann_z->Write();
      h_xy->Write("XY_denisty");
      h_yz->Write("YZ_denisty");
      h_xz->Write("XZ_denisty");
      fNew.Close();
}
