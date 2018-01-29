#include <iostream>
#include <conio.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include<iomanip>
# define reso 23000
# define co2 0.25
# define co4 0.00390625

using namespace std;

long double nx[reso], ny[reso], ro[reso], u[reso], v[reso], T[reso], p[reso],e[reso],cn[reso][5],delcn[reso][5]={0},delcn1[reso][5]={0}, c[reso][5],sedge[reso],R,F[reso][5],Fv[reso][5],G[reso][5],Gv[reso][5],mu,k,lhs[reso]={0.0},rhs[reso][5]={0.0},rrhs[reso][5]={0.0};
long double p0,T0,gamma,molwt,pr,M,cfl,cp,cv,Tin,pin,uin,vin,ein,dti,roin,dtmin=1.0e6,totime=0.0;
long double rou, rov, roeu, roev, uavg, vavg, txx[reso], tyy[reso], txy[reso], rou2, rov2, rouv, pel,pu,pv;
long double b1[reso], b2[reso], b3[reso], c1[reso], c2[reso], c3[reso], vort[reso], evort[reso];
long double dudx[reso], dudy[reso], dvdx[reso], dvdy[reso], dTdx[reso], dTdy[reso];
long double nu[reso]={0}, d2phi[reso][5], d4phi[reso][5], d24[reso][5];
long double elarea[reso], /*minidist=100.00, dist,*/ l1 , l2, l3,scx,scy,dnr,enx,eny,uav,vav,rouav,rovav,rou2av,rov2av,rouvav,rouhav,rovhav,Fav[5],Gav[5], pav;

int nbel[reso][3], nside[reso][3], rstrt, thoustep=0, totbel, bn1, bn2, bes[reso][3]={0}, timesteps;
int n1[reso], n2[reso], n3[reso],nodes, neles, neighborno[reso]={0}, neighbornode[reso][20]={0},closeneig[reso][3],ntype[reso],totbnodes;

long double length(long double x1,long double x2,long double y1,long double y2)
{long double l=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
return(l);}

long double mod(long double i)
	{ if(i<0)i*=-1; return(i);}


void tecplot()

   {long double M;

       ofstream tecfile("tec.dat");

         tecfile<<"variables ="<<"\t"<<"x,"<<"y,"<<"u,"<<"v,"<<"p, T, M, vort"<<endl;
         tecfile<<"ZONE F=FEPOINT,ET=TRIANGLE,N ="<< nodes<<" ,E = "<<neles<<endl;

       for(int i=1;i<=nodes;i++)
       {
       M=sqrt(u[i]*u[i]+v[i]*v[i])/sqrt(gamma*R*T[i]);
       tecfile<<nx[i]<<'\t'<<ny[i]<<'\t'<<u[i]<<'\t'<<v[i]<<'\t'<<(p[i]-1.0e5)<<'\t'<<T[i]<<'\t'<<M<<'\t'<<vort[i]<<endl;
       }

     for(int i=1;i<=neles;i++)

      {
         tecfile<<n1[i]<<"\t"<<n2[i]<<"\t"<<n3[i]<<endl;

      }

   }


void timestep()
	{
   	dtmin=1.0e6;

   	for(int i=1; i<=neles; i++){dti=sedge[i]/(sqrt(pow((u[n1[i]]+u[n2[i]]+u[n3[i]])/3.0,2)+pow((v[n1[i]]+v[n2[i]]+v[n3[i]])/3.0,2))+sqrt(gamma*R*(T[n1[i]]+T[n2[i]]+T[n3[i]])/3.0));
      									if(dtmin>dti)dtmin=dti;}

      dtmin*=cfl;
   }


void initialization()
	{

      R=8314.0/molwt;
      cp=R*gamma/(gamma-1.0);
      cv=cp-R;
      Tin=T0/*/(1+(gamma-1.0)*M*M*0.5)*/;
      pin=p0/*/pow((1+(gamma-1.0)*M*M*0.5),(gamma/(gamma-1.0)))*/;
      uin=sqrt(gamma*R*Tin)*M;
      vin=0.0;
      ein=cv*Tin+(pow(uin,2)+pow(vin,2))*0.5;
      roin=pin/(R*Tin);

      if(rstrt==0)for(int i=1; i<=nodes; i++)
      {
      	ro[i]=roin;
      	u[i]=uin;
      	v[i]=vin;
         p[i]=pin;

      	T[i]=Tin;
         e[i]=ein;

         cn[i][1]=ro[i];
         cn[i][2]=ro[i]*u[i];
         cn[i][3]=ro[i]*v[i];
         cn[i][4]=ro[i]*e[i];


      }

      //restart
      if (rstrt>0){
      					ifstream restartfile("restart.in");
                     for(int i=1; i<=nodes; i++)
                     {
                     	restartfile>>cn[i][1]>>cn[i][2]>>cn[i][3]>>cn[i][4]>>delcn[i][1]>>delcn[i][2]>>delcn[i][3]>>delcn[i][4];
                     	ro[i]=cn[i][1];
                        u[i]=cn[i][2]/cn[i][1];
                        v[i]=cn[i][3]/cn[i][1];
                        e[i]=cn[i][4]/cn[i][1];
                        T[i]=(e[i]-(u[i]*u[i]+v[i]*v[i])*0.5)/cv;
                        p[i]=ro[i]*R*T[i];
                     }
                     cout<<"restarted"<<endl;
                   }
   }

   void FGcal()
	{  //dudx=0; dvdy=0; dvdx=0; dudy=0; dTdx=0; dTdy=0;
   	long double pavg, uavg, vavg, rouh, rovh;

      for(int i=1; i<=neles; i++){dudx[i]=b1[i]*u[n1[i]]+b2[i]*u[n2[i]]+b3[i]*u[n3[i]];
											dvdy[i]=c1[i]*v[n1[i]]+c2[i]*v[n2[i]]+c3[i]*v[n3[i]];

                                 dvdx[i]=b1[i]*v[n1[i]]+b2[i]*v[n2[i]]+b3[i]*v[n3[i]];
                                 dudy[i]=c1[i]*u[n1[i]]+c2[i]*u[n2[i]]+c3[i]*u[n3[i]];

                                 dTdx[i]=b1[i]*T[n1[i]]+b2[i]*T[n2[i]]+b3[i]*T[n3[i]];
                                 dTdy[i]=c1[i]*T[n1[i]]+c2[i]*T[n2[i]]+c3[i]*T[n3[i]];

      									rouh=((ro[n1[i]]*e[n1[i]]+p[n1[i]])*u[n1[i]]+(ro[n2[i]]*e[n2[i]]+p[n2[i]])*u[n2[i]]+(ro[n3[i]]*e[n3[i]]+p[n3[i]])*u[n3[i]])/3.0;
                                 rovh=((ro[n1[i]]*e[n1[i]]+p[n1[i]])*v[n1[i]]+(ro[n2[i]]*e[n2[i]]+p[n2[i]])*v[n2[i]]+(ro[n3[i]]*e[n3[i]]+p[n3[i]])*v[n3[i]])/3.0;

     										rou=(ro[n1[i]]*u[n1[i]]+ro[n2[i]]*u[n2[i]]+ro[n3[i]]*u[n3[i]])/3.0;
                                 rov=(ro[n1[i]]*v[n1[i]]+ro[n2[i]]*v[n2[i]]+ro[n3[i]]*v[n3[i]])/3.0;

                                 rou2=(ro[n1[i]]*u[n1[i]]*u[n1[i]]+ro[n2[i]]*u[n2[i]]*u[n2[i]]+ro[n3[i]]*u[n3[i]]*u[n3[i]])/3.0;
                                 rov2=(ro[n1[i]]*v[n1[i]]*v[n1[i]]+ro[n2[i]]*v[n2[i]]*v[n2[i]]+ro[n3[i]]*v[n3[i]]*v[n3[i]])/3.0;
											rouv=(ro[n1[i]]*u[n1[i]]*v[n1[i]]+ro[n2[i]]*u[n2[i]]*v[n2[i]]+ro[n3[i]]*u[n3[i]]*v[n3[i]])/3.0;


                                 pavg=(p[n1[i]]+p[n2[i]]+p[n3[i]])/3.0;
                                 uavg=(u[n1[i]]+u[n2[i]]+u[n3[i]])/3.0;
                                 vavg=(v[n1[i]]+v[n2[i]]+v[n3[i]])/3.0;

                                 //sutherland
                                 mu=(11.848e-08)*sqrt(molwt)*pow((T[n1[i]]+T[n2[i]]+T[n3[i]])/3.0,0.6);

                                 txx[i]= mu*(2.0*dudx[i]-dvdy[i])*(2.0/3.0);
                                 tyy[i]= mu*(2.0*dvdy[i]-dudx[i])*(2.0/3.0);
                                 txy[i]= mu*(dvdx[i]+dudy[i]);


                                 F[i][1]=rou;									  				G[i][1]=rov;
                                 F[i][2]=rou2+pavg-txx[i];									G[i][2]=rouv-txy[i];
                                 F[i][3]=rouv-txy[i];											G[i][3]=rov2+pavg-tyy[i];
                                 F[i][4]=rouh-k*dTdx[i]-uavg*txx[i]-vavg*txy[i]; 	G[i][4]=rovh-k*dTdy[i]-vavg*tyy[i]-uavg*txy[i];

                             };
      }



void artidiss()
	{  long double nunr, nudr, c1, c4, lhsbt;

   	for(int i=1;i<=nodes;i++)
      {	nunr=0.0; nudr=0.0;

         for(int j=1; j<=neighborno[i];j++)
         {
         	nunr+=mod(p[neighbornode[i][j]]-p[i]);
         	nudr+=mod(p[neighbornode[i][j]]+p[i]);
         }


         nu[i]=nunr/nudr;
      }

      for(int k=1; k<=4; k++)
      {
      	for(int i=1;i<=nodes;i++)
         {d2phi[i][k]=0.0;

         	for(int j=1; j<=neighborno[i]; j++)
            {
             lhsbt=(lhs[neighbornode[i][j]]/dtmin+lhs[i]/dtmin)*0.5;
             d2phi[i][k]+=lhsbt*co2*max(nu[neighbornode[i][j]],nu[i])*(cn[neighbornode[i][j]][k]-cn[i][k]);
            }
         }

      }

     for(int k=1; k<=4; k++)
      {
      	for(int i=1;i<=nodes;i++)
         {d4phi[i][k]=0.0;

         	for(int j=1; j<=neighborno[i]; j++)
            {
             lhsbt=(lhs[neighbornode[i][j]]/dtmin+lhs[i]/dtmin)*0.5;
             c1=co2*max(nu[neighbornode[i][j]],nu[i]);
             if ((co4-c1)>0.0)c4=(co4-c1); else c4=0.0;
             d4phi[i][k]+=lhsbt*c4*(d2phi[neighbornode[i][j]][k]-d2phi[i][k]);

            }
         }

      }


      for(int k=1; k<=4; k++)
      {
      	for(int i=1; i<=nodes; i++)
         {
         	/*c1=co2*nu[i];
            if ((co4-c1)>0.0)c4=(co4-c1); else c4=0.0;*/
            d24[i][k]=(d2phi[i][k]-d4phi[i][k]);
         }
      }

    }

void boundcond()
	{

      for(int i=1; i<=totbnodes; i++)
      {
      	//1=supersonic inflow
         if(ntype[i]==1){  u[closeneig[i][1]]=uin;
                           v[closeneig[i][1]]=vin;
                           T[closeneig[i][1]]=Tin;
                           ro[closeneig[i][1]]=roin;
         					}

         //2=subsonic inflow
      	if(ntype[i]==2){	u[closeneig[i][1]]=uin;//sqrt(pow(u[closeneig[i][2]],2)+pow(v[closeneig[i][2]],2));
                           v[closeneig[i][1]]=vin;//0.0;
                           p[closeneig[i][1]]=pin;
                           T[closeneig[i][1]]=T0-(pow(u[closeneig[i][2]],2)+pow(v[closeneig[i][2]],2))*0.5/cp;
                           ro[closeneig[i][1]]=p[closeneig[i][1]]/(R*T[closeneig[i][1]]);
         					}

         //3=outflow
         if(ntype[i]==3){
                           u[closeneig[i][1]]=mod(u[closeneig[i][2]]);
                           v[closeneig[i][1]]=v[closeneig[i][2]];
                           p[closeneig[i][1]]=p[closeneig[i][2]];
                           T[closeneig[i][1]]=T[closeneig[i][2]];
                           ro[closeneig[i][1]]=p[closeneig[i][1]]/(R*T[closeneig[i][1]]);
         					}

         //4=symmetry
         if(ntype[i]==4){  u[closeneig[i][1]]=u[closeneig[i][2]];
         						v[closeneig[i][1]]=0.0;
                           ro[closeneig[i][1]]=ro[closeneig[i][2]];
                           p[closeneig[i][1]]=p[closeneig[i][2]];
                           T[closeneig[i][1]]=T[closeneig[i][2]];
         					}

         //5=no slip wall
         if(ntype[i]==5){  u[closeneig[i][1]]=0.0;
                          	v[closeneig[i][1]]=0.0;
                           p[closeneig[i][1]]=p[closeneig[i][2]];
                           T[closeneig[i][1]]=T[closeneig[i][2]];
                           ro[closeneig[i][1]]=p[closeneig[i][1]]/(R*T[closeneig[i][1]]);
         					}

      	e[closeneig[i][1]]=cv*T[closeneig[i][1]]+(pow(u[closeneig[i][1]],2)+pow(v[closeneig[i][1]],2))*0.5;
         //cout<<u[closeneig[i][1]]<<'\t'<<v[closeneig[i][1]]<<endl;
      }
   }

   void vorticity()
   {
   	for(int i=1; i<=neles; i++)evort[i]=b1[i]*v[n1[i]]-c1[i]*u[n1[i]]+b2[i]*v[n2[i]]-c2[i]*u[n2[i]]+b3[i]*v[n3[i]]-c3[i]*u[n3[i]];
		for(int i=1; i<=neles; i++){ vort[n1[i]]+=evort[i]/3.0;	vort[n2[i]]+=evort[i]/3.0;	vort[n3[i]]+=evort[i]/3.0;}

   }


   void results()
   {	ofstream resfile("6.57 results.dat");
     	for (int i=1; i<=nodes; i++)resfile<<i<<"\t\t"<<p[i]<<"\t\t"<<ro[i]<<"\t\t"<<T[i]<<"\t\t"<<u[i]<<"\t\t"<<v[i]<<endl;
     	//resfile<<"node\t\tP\t\tro\t\tT\t\tu\t\tv\n";
   }


int main()
	{
   int n, elno, dummy, cylnodes, innodes, outnodes, wallnodes, counter=0, z=1, boundnodes[1000][20]={0}, m=2;


   //long double rhs1, rhs2, rhs3;
   ofstream vis("visual.dat");

   ifstream gridfile("von karmann.neu");            //reading no. of nodes and elements
   gridfile>>nodes>>neles;
   cout<<nodes<<'\t'<<neles<<endl;
   vis<<nodes <<'\t'<<neles<<'\t'<<'3'<<endl<<endl;


   for(int i=1;i<=nodes;i++)                     //loop to read nodes and co-ordinates
   	{
    	gridfile>>n;
    	gridfile>>nx[n]>>ny[n];
      vis<<n<<'\t'<<nx[n]<<'\t'<<ny[n]<<'\t'<<"0.00"<<endl;

      nx[n]*=0.001;
      ny[n]*=0.001;


   	}

      ofstream writestart("element properties.dat");
      writestart<<"element no.\tn1\tn2\tn3\tb1\tb2\tb3\tc1\tc2\tc3\telement area\tshortest edge length\n";

      ofstream nenode("neighbor nodes.dat");
      nenode<<"node\tno\tneighbouring nodes\n";

      //reading boundary nodes
      ifstream cylfile("wall.dat"); 	cylfile>>dummy>>cylnodes>>dummy>>dummy;
		ifstream infile("inlet.dat"); 		infile>>dummy>>innodes>>dummy>>dummy;
      ifstream outfile("outlet.dat"); 		outfile>>dummy>>outnodes>>dummy>>dummy;
      //ifstream wallfile("wall.dat");      wallfile>>dummy>>wallnodes>>dummy>>dummy;




   for (int j=1;j<=neles;j++)                   //loop to read elements and write results
   	{
      gridfile>>elno>>dummy>>dummy>>n1[j]>>n2[j]>>n3[j];
      vis<<elno<<'\t'<<n1[j]<<'\t'<<n2[j]<<'\t'<<n3[j]<<'\t'<<'0'<<endl;



      //shape gradients calculation across elements
      elarea[j]=nx[n2[j]]*ny[n3[j]]-nx[n3[j]]*ny[n2[j]]-nx[n1[j]]*ny[n3[j]]+nx[n3[j]]*ny[n1[j]]+nx[n1[j]]*ny[n2[j]]-nx[n2[j]]*ny[n1[j]];
      elarea[j]*=0.5;

      b1[elno]=(ny[n2[j]]-ny[n3[j]])/(2*elarea[j]);	c1[elno]=(nx[n3[j]]-nx[n2[j]])/(2*elarea[j]);
      b2[elno]=(ny[n3[j]]-ny[n1[j]])/(2*elarea[j]);	c2[elno]=(nx[n1[j]]-nx[n3[j]])/(2*elarea[j]);
      b3[elno]=(ny[n1[j]]-ny[n2[j]])/(2*elarea[j]);	c3[elno]=(nx[n2[j]]-nx[n1[j]])/(2*elarea[j]);



      //finding out the shortest edge of an element
      l1=length(nx[n1[j]],nx[n2[j]],ny[n1[j]],ny[n2[j]]);
      l2=length(nx[n2[j]],nx[n3[j]],ny[n2[j]],ny[n3[j]]);
      l3=length(nx[n3[j]],nx[n1[j]],ny[n3[j]],ny[n1[j]]);

      sedge[elno]=l1; if(l2<sedge[elno])sedge[elno]=l2; if(l3<sedge[elno])sedge[elno]=l3;

      writestart<<elno<<'\t'<<n1[j]<<'\t'<<n2[j]<<'\t'<<n3[j]<<'\t'<<b1[elno]<<'\t'<<b2[elno]<<'\t'<<b3[elno]<<'\t'<<c1[elno]<<'\t'<<c2[elno]<<'\t'<<c3[elno]<<'\t'<<elarea[elno]<<'\t'<<sedge[elno]<<endl;

      lhs[n1[j]]+=elarea[j]/3.0;
      lhs[n2[j]]+=elarea[j]/3.0;
      lhs[n3[j]]+=elarea[j]/3.0;

      //finding out neighboring nodes to every node
      for(int i=1;i<=nodes;i++)	{

      											if(i==n1[j]){	while(neighbornode[i][z]!=0){if(neighbornode[i][z]==n2[j])counter++;z++;};     			//checking if n2 is already written.
                                       																																	//all elements of neighbornode are initalized to 0.
                                                                                                                  								//all elements checked if =n2. if counter=0, then first unfilled element =n2.
                                                                  if(counter==0)neighbornode[i][z]=n2[j]; z=1; counter=0;


                                                   					while(neighbornode[i][z]!=0){if(neighbornode[i][z]==n3[j])counter++;z++;};    //checzing if n3 is already written
                                                                  																										//same algo for n3.

                                                                  if(counter==0)neighbornode[i][z]=n3[j]; z=1; counter=0;
                                       								};

                                       if(i==n2[j]){	while(neighbornode[i][z]!=0){if(neighbornode[i][z]==n1[j])counter++;z++;};


                                                                  if(counter==0)neighbornode[i][z]=n1[j];z=1; counter=0;


                                                                  while(neighbornode[i][z]!=0){if(neighbornode[i][z]==n3[j])counter++;z++;};


                                                                  if(counter==0)neighbornode[i][z]=n3[j];z=1; counter=0;
                                       								};

                                       if(i==n3[j]){	while(neighbornode[i][z]!=0){if(neighbornode[i][z]==n2[j])counter++;z++;};


                                                                  if(counter==0)neighbornode[i][z]=n2[j];z=1; counter=0;


                                                   					while(neighbornode[i][z]!=0){if(neighbornode[i][z]==n1[j])counter++;z++;};


                                                                  if(counter==0)neighbornode[i][z]=n1[j]; counter=0;z=1;
                                                 					};

                                 	while(neighbornode[i][++neighborno[i]]!=0);
                                    neighborno[i]--; //cout<<i<<' '<<neighborno[i]<<'\t';
      									};
      }


       for(int i=1;i<=nodes;i++){for(int j=1;j<=10;j++){if(j==1)nenode<<i<<'\t'<<neighborno[i]; if(neighbornode[i][j]!=0)nenode<<'\t'<<neighbornode[i][j];}nenode<<endl;}

      //reading all boundary nodes from different files to one array
      for(int i=1;i<=cylnodes;i++){cylfile>>boundnodes[i][1];ntype[i]=5;}
      for(int i=cylnodes+1;i<=cylnodes+innodes;i++){infile>>boundnodes[i][1];ntype[i]=1;}
      for(int i=cylnodes+innodes+1;i<=cylnodes+innodes+outnodes;i++){outfile>>boundnodes[i][1];ntype[i]=3;}


      totbnodes=cylnodes+innodes+outnodes;

      //finding out the neighbor nodes of boundary nodes
      counter=0;
      for(int i=1;i<=totbnodes;i++)
      {for(int j=1;j<=nodes;j++)if(boundnodes[i][1]==neighbornode[j][1])
      																						for(int k=1; k<=neighborno[j]; k++)
                                                                        												{for(int l=1;l<=totbnodes;l++)if(neighbornode[j][k]!=boundnodes[l][1])counter++;if(counter==totbnodes)boundnodes[i][m++]=neighbornode[j][k];counter=0;}m=2;}

   	ifstream bcfile("bc.in");
      //bcfile>>dummy;
      for(int i=1;i<=totbnodes;i++)bcfile>>closeneig[i][1]>>closeneig[i][2]>>ntype[i]>>nbel[i][1]>>nside[i][1]>>nbel[i][2]>>nside[i][2];

      //reading all boundary elements to a file

      counter=0; int l=1;
      for(int i=1;i<=totbnodes;i++){for(int j=1;j<=l;j++)if(bes[j][1]==nbel[i][1]&&bes[j][2]==nside[i][1])counter++; if(counter==0){ bes[l][1]=nbel[i][1]; bes[l++][2]=nside[i][1];}counter=0;}
      for(int i=1;i<=totbnodes;i++){for(int j=1;j<=l;j++)if(bes[j][1]==nbel[i][2]&&bes[j][2]==nside[i][2])counter++; if(counter==0){ bes[l][1]=nbel[i][2]; bes[l++][2]=nside[i][2];}counter=0;}

     totbel=--l;

      ifstream flowfile("flow.in");
      flowfile>>p0>>T0>>gamma>>molwt>>pr>>M>>mu>>k>>cfl>>timesteps>>rstrt;

      initialization();
     	boundcond();
     	artidiss();
      timestep();
      FGcal();


 		//time step iterations
for(int tstep=1; tstep<=timesteps; tstep++)
{
      for (int i=1; i<=neles; i++)
      {
      	for (int j=1; j<=4; j++)
         {
         	//evaluating RHS
      		rhs[n1[i]][j]+=elarea[i]*(F[i][j]*b1[i]+G[i][j]*c1[i]);
            rhs[n2[i]][j]+=elarea[i]*(F[i][j]*b2[i]+G[i][j]*c2[i]);
            rhs[n3[i]][j]+=elarea[i]*(F[i][j]*b3[i]+G[i][j]*c3[i]);
            }

            //boundary treatment
         for(int m=1; m<=totbel;m++)if(i==bes[m][1]){
            															if(bes[m][2]==1){bn1=n1[i];bn2=n2[i];}
                                                      	if(bes[m][2]==2){bn1=n2[i];bn2=n3[i];}
                                                      	if(bes[m][2]==3){bn1=n3[i];bn2=n1[i];}

            															scx=ny[bn2]-ny[bn1];
            															scy=nx[bn1]-nx[bn2];
            															dnr=sqrt(scx*scx+scy*scy);

                                                      	enx=scx/dnr;
                                                      	eny=scy/dnr;

                                                      	uav=(u[bn1]+u[bn2])*0.5;												vav=(v[bn1]+v[bn2])*0.5;
	                                                      rouav=(ro[bn1]*u[bn1]+ro[bn2]*u[bn2])*0.5;						rovav=(ro[bn1]*v[bn1]+ro[bn2]*v[bn2])*0.5;
   	                                                   rou2av=(ro[bn1]*u[bn1]*u[bn1]+ro[bn2]*u[bn2]*u[bn2])*0.5;	rov2av=(ro[bn1]*v[bn1]*v[bn1]+ro[bn2]*v[bn2]*v[bn2])*0.5;
      	                                                rouvav=(ro[bn1]*u[bn1]*v[bn1]+ro[bn2]*u[bn2]*v[bn2])*0.5;   pav=(p[bn1]+p[bn2])*0.5;
         	                                             rouhav=((ro[bn1]*e[bn1]+p[bn1])*u[bn1]+(ro[bn2]*e[bn2]+p[bn2])*u[bn2])*0.5;
            	                                          rovhav=((ro[bn1]*e[bn1]+p[bn1])*v[bn1]+(ro[bn2]*e[bn2]+p[bn2])*v[bn2])*0.5;

                                                         Fav[1]=rouav;						Gav[1]=rovav;
                  	                                    Fav[2]=rou2av+pav-txx[i];		Gav[2]=rouvav-txy[i];
                     	                                 Fav[3]=rouvav-txy[i];			Gav[3]=rov2av+pav-tyy[i];
                        	                              Fav[4]=rouhav-k*dTdx[i]-uav*txx[i]-vav*txy[i];
                           	                           Gav[4]=rovhav-k*dTdy[i]-vav*tyy[i]-uav*txy[i];

                                                         	for(int j=1; j<=4; j++)
                                                         	{
            																	rhs[bn1][j]-=dnr*0.5*(Fav[j]*enx+Gav[j]*eny);
                                                               rhs[bn2][j]-=dnr*0.5*(Fav[j]*enx+Gav[j]*eny);
                                                         	}
            														  }



      }


      if(tstep==1&& rstrt==0) for(int i=1; i<=nodes; i++)for(int j=1; j<=4; j++){delcn[i][j]=dtmin*(rhs[i][j]+d24[i][j])/(3.0*lhs[i]);rhs[i][j]=0.0;}


      // convergence iteration
     if(tstep>1||rstrt==1)for(int iter=1; iter<=10; iter++)
     {
     		//evaluating RRHS
     		for(int i=1; i<=neles; i++)for(int j=1; j<=4; j++)
     		{	rrhs[n1[i]][j]+=(2*delcn[n1[i]][j]-delcn[n2[i]][j]-delcn[n3[i]][j])*elarea[i]/12.0;
     			rrhs[n2[i]][j]+=(-delcn[n1[i]][j]+2*delcn[n2[i]][j]-delcn[n3[i]][j])*elarea[i]/12.0;
        	   rrhs[n3[i]][j]+=(-delcn[n1[i]][j]-delcn[n2[i]][j]+2*delcn[n3[i]][j])*elarea[i]/12.0;
     		}

     		for(int i=1; i<=nodes; i++)for(int j=1; j<=4; j++){
         																	delcn1[i][j]=(dtmin*(rhs[i][j]+d24[i][j])+rrhs[i][j])/(3.0*lhs[i]);
         																	rrhs[i][j]=0.0;
																				}

         //in-iteration error calc
         //cout<<tstep<<'\t'<<iter<<'\t'<<(delcn1[56][1]-delcn[56][1])<<'\t'<<(delcn1[56][2]-delcn[56][2])<<'\t'<<(delcn1[56][3]-delcn[56][3])<<'\t'<<(delcn1[56][4]-delcn[56][4])<<endl;

      	for(int i=1; i<=nodes; i++)for (int j=1; j<=4; j++){delcn[i][j]=delcn1[i][j]; delcn1[i][j]=0.0;}


     }
     totime+=dtmin;
     //cout<<tstep<<'\t'<<totime<<'\t'<<delcn[11470][1]<<'\t'<< delcn[11470][2]<<'\t'<< delcn[11470][3]<<'\t'<< delcn[11470][4]<<endl;




      for(int i=1; i<=nodes; i++)
      {
         u[i]=(delcn[i][2]+ro[i]*u[i])/(ro[i]+delcn[i][1]);
      	v[i]=(delcn[i][3]+ro[i]*v[i])/(ro[i]+delcn[i][1]);
      	e[i]=(delcn[i][4]+ro[i]*e[i])/(ro[i]+delcn[i][1]);
         T[i]=(e[i]-(pow(u[i],2)+pow(v[i],2))*0.5)/cv; //if (T[i]<=0)T[i]=1;
         ro[i]+=delcn[i][1];
         p[i]=ro[i]*R*T[i];

      	cn[i][1]=ro[i];
         cn[i][2]=ro[i]*u[i];
			cn[i][3]=ro[i]*v[i];
         cn[i][4]=ro[i]*e[i];
      }

      for(int i=1; i<=nodes; i++)for(int j=1; j<=4; j++)rhs[i][j]=0.0;
	  
	  
      //tecplot();
      FGcal();
      artidiss();
      boundcond();
      timestep();
	  
	  //cout<<tstep<<'\t'<<totime<<'\t'<<u[25]<<endl;
	  
      thoustep++;
      if(thoustep==1000){
      						ofstream refile("restart.in");
      						for(int i=1; i<=nodes; i++)refile<<cn[i][1]<<"\t\t"<<cn[i][2]<<"\t\t"<<cn[i][3]<<"\t\t"<<cn[i][4]<<"\t\t"<<delcn[i][1]<<"\t\t"<<delcn[i][2]<<"\t\t"<<delcn[i][3]<<"\t\t"<<delcn[i][4]<<endl;
                        vorticity(); tecplot(); results(); thoustep=0;
                        cout<<tstep<<'\t'<<totime<<'\t'<<delcn[11470][1]<<'\t'<< delcn[11470][2]<<'\t'<< delcn[11470][3]<<'\t'<< delcn[11470][4]<<endl;
                        }

}

		ofstream refile("restart.in");
      for(int i=1; i<=nodes; i++)refile<<setprecision(8)<<cn[i][1]<<"\t\t"<<setprecision(8)<<cn[i][2]<<"\t\t"<<setprecision(8)<<cn[i][3]<<"\t\t"<<setprecision(8)<<cn[i][4]<<"\t\t"<<setprecision(8)<<delcn[i][1]<<"\t\t"<<setprecision(8)<<delcn[i][2]<<"\t\t"<<setprecision(8)<<delcn[i][3]<<"\t\t"<<setprecision(8)<<delcn[i][4]<<endl;

		boundcond();
      tecplot();
      results();

     getch();
   }
