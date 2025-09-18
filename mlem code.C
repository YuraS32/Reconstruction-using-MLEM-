#include <RtypesCore.h>
#include <algorithm>
#define medicalProject_cxx
#include "medicalProject.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <math.h>
#include <TGraph2D.h>
#include <TH3D.h>

struct VoxelSeg
{
    int i,j,k;
    double len;
};

void medicalProject::Loop()
{
    if (fChain == 0) return;

    TFile fNew("output_8iter.root","RECREATE");

    double Xmin=-20, Xmax=20;
    double Ymin=-20, Ymax=20;
    double Zmin=-20, Zmax=20;
    double dx = 0.2, dy = 0.2, dz = 0.5;
    int Nx = (Xmax-Xmin)/dx;
    int Ny = (Ymax-Ymin)/dy;
    int Nz = (Zmax-Zmin)/dz;

    std::vector<float> f(Nx*Ny*Nz, 1.0f);
    std::vector<float> delta(Nx*Ny*Nz, 0.0f);
    std::vector<float> sumP(Nx*Ny*Nz, 0.0f);

    auto idx = [&](int i,int j,int k){ return (i*Ny + j)*Nz + k; };

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    for(Long64_t jentry=0;jentry<nentries;jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if(ientry<0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        double dx_hit = X_hit2-X_hit1;
        double dy_hit = Y_hit2-Y_hit1;
        double dz_hit = Z_hit2-Z_hit1;
        double L = sqrt(dx_hit*dx_hit + dy_hit*dy_hit + dz_hit*dz_hit);
        if(L<=0.0f) continue;

        std::vector<double> v;
        if(fabs(dx_hit)>1e-12)
            for(double x=Xmin+dx;x<Xmax;x+=dx)
            {
                double t=(x-X_hit1)/dx_hit;
                if(t>1e-9 && t<1.0-1e-9) v.push_back(t);
            }

        if(fabs(dy_hit)>1e-12)
            for(double y=Ymin+dy;y<Ymax;y+=dy)
            {
                double t=(y-Y_hit1)/dy_hit;
                if(t>1e-9 && t<1.0-1e-9) v.push_back(t);
            }

        if(fabs(dz_hit)>1e-12)
            for(double z=Zmin+dz;z<Zmax;z+=dz)
            {
                double t=(z-Z_hit1)/dz_hit;
                if(t>1e-9 && t<1.0-1e-9) v.push_back(t);
            }

        v.push_back(0.0);
        v.push_back(1.0);
        std::sort(v.begin(),v.end());
        v.erase(std::unique(v.begin(),v.end()),v.end());

        for(size_t i=0;i+1<v.size();i++)
        {
            double t_mid=0.5*(v[i]+v[i+1]);
            double x_v=X_hit1+dx_hit*t_mid;
            double y_v=Y_hit1+dy_hit*t_mid;
            double z_v=Z_hit1+dz_hit*t_mid;
            int i_v=floor((x_v-Xmin)/dx);
            int j_v=floor((y_v-Ymin)/dy);
            int k_v=floor((z_v-Zmin)/dz);
            if(i_v>=0 && i_v<Nx && j_v>=0 && j_v<Ny && k_v>=0 && k_v<Nz)
                sumP[idx(i_v,j_v,k_v)]+=(v[i+1]-v[i])*L;
        }
    }

    f=sumP;

    int m=6;

    for(int it=0; it<m; it++)
    {
        std::fill(delta.begin(),delta.end(),0.0f);
        int flag=0;

        for(Long64_t jentry=0;jentry<nentries;jentry++)
        {
            flag++;
            if(flag%1000==0) std::cout<<"iteration: "<<it+1<<", event: "<<flag<<"\n";

            Long64_t ientry=LoadTree(jentry);
            if(ientry<0) break;
            nb=fChain->GetEntry(jentry);
            nbytes+=nb;

            double dx_hit=X_hit2-X_hit1;
            double dy_hit=Y_hit2-Y_hit1;
            double dz_hit=Z_hit2-Z_hit1;
            double L=sqrt(dx_hit*dx_hit + dy_hit*dy_hit + dz_hit*dz_hit);

            if(L<=0.0f) continue;

            std::vector<double> v;
            if(fabs(dx_hit)>1e-12)
                for(double x=Xmin+dx;x<Xmax;x+=dx)
                {
                    double t=(x-X_hit1)/dx_hit;
                    if(t>1e-9 && t<1.0-1e-9) v.push_back(t);
                }

            if(fabs(dy_hit)>1e-12)
                for(double y=Ymin+dy;y<Ymax;y+=dy)
                {
                    double t=(y-Y_hit1)/dy_hit;
                    if(t>1e-9 && t<1.0-1e-9) v.push_back(t);
                }

            if(fabs(dz_hit)>1e-12)
                for(double z=Zmin+dz;z<Zmax;z+=dz)
                {
                    double t=(z-Z_hit1)/dz_hit;
                    if(t>1e-9 && t<1.0-1e-9) v.push_back(t);
                }

            v.push_back(0.0);
            v.push_back(1.0);
            std::sort(v.begin(),v.end());
            v.erase(std::unique(v.begin(),v.end()),v.end());

            double s_j=0.0;
            struct Tmp{int i,j,k; double len;};
            std::vector<Tmp> voxels;

            for(size_t mseg=0; mseg+1<v.size(); mseg++)
            {
                double t_mid=0.5*(v[mseg]+v[mseg+1]);
                double x_v=X_hit1+dx_hit*t_mid;
                double y_v=Y_hit1+dy_hit*t_mid;
                double z_v=Z_hit1+dz_hit*t_mid;
                int i_v=floor((x_v-Xmin)/dx);
                int j_v=floor((y_v-Ymin)/dy);
                int k_v=floor((z_v-Zmin)/dz);
                if(i_v>=0 && i_v<Nx && j_v>=0 && j_v<Ny && k_v>=0 && k_v<Nz)
                {
                    double len=(v[mseg+1]-v[mseg])*L;
                    voxels.push_back({i_v,j_v,k_v,len});
                    s_j+=len*f[idx(i_v,j_v,k_v)];
                }
            }

            if(s_j<1e-12) continue;

            for(auto &vt:voxels)
                delta[idx(vt.i,vt.j,vt.k)]+=vt.len*f[idx(vt.i,vt.j,vt.k)]/s_j;
            
        }

        for(int i=0;i<Nx;i++)
            for(int j=0;j<Ny;j++)
                for(int k=0;k<Nz;k++)
                    if(sumP[idx(i,j,k)]>1e-12)
                        f[idx(i,j,k)] *= std::max(delta[idx(i,j,k)]/sumP[idx(i,j,k)], 1e-3f);


        TH3D *h3=new TH3D("activity","MLEM activity",Nx,Xmin,Xmax,Ny,Ymin,Ymax,Nz,Zmin,Zmax);
    for(int i=0;i<Nx;i++)
        for(int j=0;j<Ny;j++)
            for(int k=0;k<Nz;k++)
                h3->SetBinContent(i+1,j+1,k+1,f[idx(i,j,k)]);

    h3->Write();
    h3->GetXaxis()->SetRangeUser(-18,18);
    h3->GetYaxis()->SetRangeUser(-18,18);
    h3->GetZaxis()->SetRangeUser(-18,18);

int nBinsZ = h3->GetNbinsZ();
int nBinsX = h3->GetNbinsX();
int nBinsY = h3->GetNbinsY();

double xlow = h3->GetXaxis()->GetBinLowEdge(1);
double xup  = h3->GetXaxis()->GetBinUpEdge(nBinsX);
double ylow = h3->GetYaxis()->GetBinLowEdge(1);
double yup  = h3->GetYaxis()->GetBinUpEdge(nBinsY);

int step = std::max(1, (int)(5.0 / dz + 0.5));

for (int k = 1; k <= nBinsZ; k += step) {
    int kmax = std::min(k + step - 1, nBinsZ);

    double zlow  = h3->GetZaxis()->GetBinLowEdge(k);
    double zhigh = h3->GetZaxis()->GetBinUpEdge(kmax);

    TString name  = Form("projYX_z%d_%d", k, kmax); 
    TString title = Form("Projection YX, Z = %.2f..%.2f cm", zlow, zhigh);

    TH2D *proj = new TH2D(name.Data(), title.Data(),
                          nBinsY, ylow, yup,   
                          nBinsX, xlow, xup);  

    for (int ix = 1; ix <= nBinsX; ++ix) {
        for (int iy = 1; iy <= nBinsY; ++iy) {
            double sum = 0.0;
            for (int iz = k; iz <= kmax; ++iz) {
                sum += h3->GetBinContent(ix, iy, iz);
            }
            proj->SetBinContent(iy, ix, sum);
        }
    }

    proj->GetXaxis()->SetRangeUser(-18, 18);
    proj->GetYaxis()->SetRangeUser(-18, 18);

    proj->Write();
    delete proj; 
}


    TH2D *projXY=(TH2D*)h3->Project3D("xy");
    projXY->Write();
    TH2D *projYZ=(TH2D*)h3->Project3D("yz");
    projYZ->Write();
    TH2D *projXZ=(TH2D*)h3->Project3D("xz");
    projXZ->Write();            
                           
}

    

    fNew.Close();
}


