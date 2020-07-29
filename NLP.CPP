//**************************************************************************
//	INCLUDED HEADER FILES
//**************************************************************************
#include<fstream.h>
#include<conio.h>
#include<stdio.h>
#include<ctype.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<iomanip.h>
#include<bios.h>
#include<dos.h>
//***************************************************************************
//	PHYSICAL CONSTANTS
//***************************************************************************
double k= 9.0e9;
double pi= M_PI;
double eps0= 8.854e-12;
double e= 1.6e-19;
double N_a=6.022e23;
//**************************************************************************
//	GLOBAL CLASSES
//**************************************************************************
struct prefix
{	char ch[4];
	int exp;
};

struct signature
{	int sign[18],total,n;
	//F,q,E,p,phi,r,lcd,scd,vcd,S,tau,angle,number,l,radius/edge,
	//position,shape,volume
	//M,a,re,z,density
} qsign={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

struct Case
{	int n;
	char str[30];
};

struct words
{	char keyword[30];
};

class quantity
{	double magnitude;
	char unit[11];
	public:
	void getdata(double,char[]);
	void getdata(double);
	void getdata(char[]);
	double retmag()
	{  if( magnitude<0 )
		return -magnitude;
		return magnitude;
	}
	char* retunit()
	{  return unit; }
	void print_quantity()
	{	cout<<" "<<magnitude<<" "<<unit; }
} data;
//**************************************************************************
// LOGICAL VARIABLES
//**************************************************************************
double F[5],q[5],r[5],E[5],S[5];
double &F5=F[4],&F1=F[0],&F2=F[1],&F3=F[2],&F4=F[3];
double &q5=q[4],&q1=q[0],&q2=q[1],&q3=q[2],&q4=q[3];
double &r5=r[4],&r1=r[0],&r2=r[1],&r3=r[2],&r4=r[3];
double &E5=E[4],&E1=E[0],&E2=E[1],&E3=E[2],&E4=E[3];
double &S5=S[4],&S1=S[0],&S2=S[1],&S3=S[2],&S4=S[3];
double lamda,sigma,rho,l,re,p,n,flux,torque,angle,vol;
double z,M,a,density;
//**************************************************************************
//	GLOBAL VARIABLES
//**************************************************************************
char tmp[1001],input[100][31],cstring[30],shape[10]="null";
char position[20]="null"; int itf=0,tw,to_find[5]={-2,-2,-2,-2,-2};
char d=196,sq=253,cpi=227,ceps=238,cs=229,ctau=226,ctheta=233,cphi=237;
char cdell=127,carr=26,sqr=251,line=179,pack[30]="null",mode;
//**************************************************************************
//	FUNCTION PROTOTYPES
//**************************************************************************
void getwords();
void convert_to_SI(quantity&);
void data_initializer_phy(char[]);
void data_initializer_chem(char[]);
void data_extracter_phy();
void data_extracter_chem();
void shape_identifier();
void case_selecter_phy(int);
int case_selecter_phy(char[]);
void case_selecter_chem(int);
int case_selecter_chem(char[]);
void coloumbs_law(char[]);
void line_charge(char[]);
void plane_charge(char[]);
void spherical_shell(char[]);
void surface_area();
void dipole_torque(char[]);
void dipole_field(char[]);
void gauss_law(char[]);
void point_field(char[]);
void calc_flux(char[]);
void calc_dipole_moment(char[]);
void quantisation_of_charge(char[]);
void calc_linear_density(char[]);
void calc_surface_density(char[]);
void calc_volume_density(char[]);
void force_field(char[]);
void calc_volume();
void cell_density(char[]);
void pack_identifier();
void radius_edge(char[]);
int spell_check(char[],char[]);
void aim_phy();
int aim_phy(char[]);
void aim_phy(char[],int&,int);
void aim_chem();
int aim_chem(char[]);
void aim_chem(char[],int&,int);
void remove_punctuation(char[]);
void check_dipole_length();
void function_selecter_phy();
void function_selecter_chem();
int sign_matcher_full_phy(signature);
int signmp(signature);
int sign_matcher_full_chem(signature);
int sign_matcher_inc_chem(signature);
void function_caller_phy(int);
void function_caller_chem(int);
void auto_mode_phy();
void auto_mode_chem();
void initializer();
void variable_printer();
void summary();
void help();
void scr();
void jnv();
void name();
void exitdesign();
void tone();
void about_matrix();
void operation();
void mat_sum();
void mat_sub();
void mat_multi();
void mat_tran();
void mat_inv();
void smcursor();
void pccursor();
void matcursor();
void inmatcursor();
void smmenu1();
void smmenu2();
void smmenu3();
void smmenu4();
void smmenu5();
void pcmenu1();
void pcmenu2();
void pcmenu3();
void pcmenu4();
void matmenu1();
void matmenu2();
void matmenu3();
void matmenu4();
void inmatmenu1();
void inmatmenu2();
void inmatmenu3();
void inmatmenu4();
void inmatmenu5();
void inmatmenu6();
//////////////////////
int smpos,ppos,cpos,matpos,inmatpos;
//**************************************************************************
//	FUNCTION NAME : initializer
//	DETAILS		  : IT INITIALIZES THE ALL GLOBAL VARIABLES
//***************************************************************************
void initializer()
{	for(int i=0;i<18;++i)
	qsign.sign[i]=0;
	for(i=0;i<5;++i)
	{	to_find[i]=0;
		r[i]=q[i]=S[i]=F[i]=E[i]=12321;
	}
	lamda=sigma=rho=l=re=p=n=flux=torque=angle=vol=M=density=a=12321;
	itf=0;
	strcpy(tmp,"");
}
//***************************************************************************
//	FUNCTION NAME : variable_printer
//	DETAILS		  : IT PRINTS THE VALUES OF THE VARIABLES
//***************************************************************************
void variable_printer()
{  int i;
	cout<<"\n\n The values of the given quantities after converting into SI:\n\n";
	for(i=0;i<5;++i)
	{	if( E[i]!=12321 )
		  cout<<"E"<<i+1<<"="<<E[i]<<" N_C(-1)  ";
	}
	for(i=0;i<5;++i)
	{	if( q[i]!=12321 )
			cout<<"q"<<i+1<<"="<<q[i]<<" C  ";
	}
	for(i=0;i<5;++i)
	{	if( r[i]!=12321 )
			cout<<"r"<<i+1<<"="<<r[i]<<" m  ";
	}
	for(i=0;i<5;++i)
	{	if( S[i]!=12321 )
			cout<<"S"<<i+1<<"="<<S[i]<<" m(2)  ";
	}
	for(i=0;i<5;++i)
	{	if( F[i]!=12321 )
			cout<<"F"<<i+1<<"="<<F[i]<<" N  ";
	}
	if(lamda!=12321) cout<<"lamda="<<lamda<<" C_m(-1)  ";
	if(sigma!=12321) cout<<"sigma="<<sigma<<" C_m(-2)  ";
	if(rho!=12321) cout<<"rho="<<rho<<" C_m(-3)  ";
	if(p!=12321) cout<<"p="<<p<<" C_m  ";
	if(flux!=12321) cout<<cphi<<" = "<<flux<<" V_m  ";
	if(torque!=12321) cout<<ctau<<"="<<torque<<" N_m  ";
	if(angle!=12321) cout<<ctheta<<"="<<angle<<" rad  ";
	if(n!=12321) cout<<"n="<<n<<"   ";
	if(l!=12321) cout<<"l="<<l<<" m  ";
	if(vol!=12321) cout<<"V="<<vol<<" m(3)  ";
	if(re!=12321 && !strcmp(shape,"circle") ) cout<<"radius="<<re<<" m (circle) ";
	if(re!=12321 && !strcmp(shape,"sphere") ) cout<<"radius="<<re<<" m (sphere)  ";
	if(re!=12321 && !strcmp(shape,"square") ) cout<<"edge="<<re<<" m (square)  ";
	if(re!=12321 && !strcmp(shape,"cube") ) cout<<"edge="<<re<<" m (cube)  ";
	if(M!=12321) cout<<"M="<<M<<" g_mol(-1)  ";
	if(a!=12321) cout<<"a="<<a<<" m  ";
	if(density!=12321) cout<<"d="<<density<<" kg_m(-3)  ";
	getch();
}

//**************************************************************************
//	CLASS NAME    : quantity
//	FUNCTION NAME : getdata
//	DETAILS       : IT TAKES DATA FOR DATA MEMBERS
//			OF CLASS QUANTITY
//**************************************************************************
void quantity::getdata(double a,char b[])
{	magnitude=a;
	strcpy(unit,b);
}

void quantity::getdata(double a)
{	magnitude=a;
}
void quantity::getdata(char b[])
{	strcpy(unit,b);
}
//***************************************************************************
//	FUNCTION NAME : convert_to_SI
//	DETAILS       : IT CONVERTS THE UNIT AND MAGNITUDE
//			OF A QUANTITY INTO SI SYSTEM
//***************************************************************************
void convert_to_SI(quantity& A)
{	ifstream fin("B_PREFIX.dat",ios::in|ios::binary);
	int i;
	if(!fin)
	{	cout<<" B_PREFIX.dat File not found";
		getch();
	}
	fin.seekg(0);
	prefix in;
	char cmp[5]; strcpy(cmp,A.retunit());
	while(fin)
	{	if( strcmp(cmp,"deg")==0 )
		{	A.getdata(A.retmag()*pi/180,"rad");
			break;
		}
		else if( cmp[0]=='d' && cmp[1]=='a')
		{  for(i=0;cmp[i]!='\0';++i)
			cmp[i]=cmp[i+2];
			A.getdata(A.retmag()*10,cmp);
			break;
		}
		else if( strcmp(cmp,"m")==0||strcmp(cmp,"m(2)")==0||strcmp(cmp,"m(3)")==0 )
			break;
		else if( strcmp(cmp,"g_cm(-3)")==0 )
		{  A.getdata(A.retmag()*1e3,"kg_m(-3)");
			break;
		}
		else if( strcmp(cmp,"u")==0 )
			A.getdata("g_mol");
		fin.read((char*)&in,sizeof(in));
		if(in.ch[0]==cmp[0])
		{  for(i=0;cmp[i]!='\0';i++)
				cmp[i]=cmp[i+1];
			A.getdata(A.retmag()*pow(10,in.exp),cmp);
			break;
		}

	}
	fin.close();
}
//***************************************************************************
//	FUNCTION NAME : case_selecter
//	DETAILS 		  : IT SELECTS THE REQUIRED CASE
//					  	 NEEDED TO SOLVE THE PROBLEM
//***************************************************************************
void case_selecter_phy(int match)
{  Case cselect;
	ifstream fin("B_CASE_SELECT.dat",ios::in|ios::binary);
	if(!fin)
	{	cout<<"\n B_CASE_SELECT.dat not found";
		getch();
	}
	fin.seekg(0);
	while(fin)
	{	fin.read( (char*)&cselect,sizeof(cselect) );
		if(match==cselect.n)
		strcpy(cstring,cselect.str);
	}
	fin.close();
}

int case_selecter_phy(char match[])
{  Case cselect;
	ifstream fin("B_CASE_SELECT.dat",ios::in|ios::binary);
	if(!fin)
	{	cout<<"\n B_CASE_SELECT not found";
		getch();
	}
	fin.seekg(0);
	while(fin)
	{	fin.read( (char*)&cselect,sizeof(cselect) );
		if( strcmp(match,cselect.str)==0 )
		{  fin.close();
			return cselect.n;
		}
	}
	fin.close();
	return -1;
}

void case_selecter_chem(int match)
{  Case cselect;
	ifstream fin("BC_CASE_SELECT.dat",ios::in|ios::binary);
	if(!fin)
	{	cout<<"\n BC_CASE_SELECT.dat not found";
		getch();
	}
	fin.seekg(0);
	while(fin)
	{	fin.read( (char*)&cselect,sizeof(cselect) );
		if(match==cselect.n)
		strcpy(cstring,cselect.str);
	}
	fin.close();
}

int case_selecter_chem(char match[])
{  Case cselect;
	ifstream fin("BC_CASE_SELECT.dat",ios::in|ios::binary);
	if(!fin)
	{	cout<<"\n BC_CASE_SELECT not found";
		getch();
	}
	fin.seekg(0);
	while(fin)
	{	fin.read( (char*)&cselect,sizeof(cselect) );
		if( strcmp(match,cselect.str)==0 )
		{  fin.close();
			return cselect.n;
		}
	}
	fin.close();
	return -1;
}

//***************************************************************************
// FUNCTION NAME : spell_check
//	DETAILS		  : IT COMPARES THE SPELLINGS OF TWO STRINGS
//***************************************************************************
int spell_check(char A[], char B[])
{  if( strlen(A)==strlen(B) )
	{  strupr(A); strupr(B);
		if( strcmp(A,B)==0 )
		return 1;
		else return 0;
	}
	else return 0;
}

//###########################################################################
//#################### FORMULAE FUNCTIONS ###################################
//***************************************************************************
//	FUNCTION NAME : coloumbs_law
//	DETAILS       : IT TAKES INPUT OF WHAT TO FIND
//		      	AND RETURNS IT USING COLOUMB'S LAW
//***************************************************************************
void coloumbs_law(char find[])
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"force")==0)
	{  gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";
		getch();
		 cout<<"\n\n\t\t\t    ("<<k<<") x ("<<q1<<") x ("<<q2<<')';
		 cout<<"\n\t\t\tF= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		 cout<<"\n\t\t\t\t     ("<<r1<<")"<<sq;  getch();
		F1= (k*q1*q2)/(r1*r1);
		cout<<"\n The force is:"<<F1<<" N";  getch();
	}

	else if( strcmp(find,"distance")==0)
	{  cout<<"  "<<(char)218<<d<<"\t     "<<d<<(char)191<<(char)171;
		cout<<"\n\t\t\t  "<<line<<"k x q1 x q2"<<line;
		cout<<"\n\t\t\tr="<<line<<""<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t  "<<line<<"     F     "<<line;
		cout<<"\n\t\t\t  "<<(char)192<<d<<"\t     "<<d<<(char)217;
		getch();
		cout<<"\n\n\t\t\t  "<<(char)218<<d<<"\t\t\t     "<<d<<(char)191<<(char)171;
		cout<<"\n\t\t\t  "<<line<<'('<<k<<") x ("<<q1<<") x ("<<q2<<')'<<line;
		cout<<"\n\t\t\tr="<<line<<""<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t  "<<line<<"\t     ("<<F1<<")          "<<line;
		cout<<"\n\t\t\t  "<<(char)192<<d<<"\t\t\t     "<<d<<(char)217;
		getch();
		r1= sqrt( k*q1*q2/F1);
		cout<<"\n The distance is:"<<r1<<" m"; getch();
	}
	else if( strcmp(find,"charge")==0)
	{  cout<<"     F x r"<<sq;
		cout<<"\n\t\t\tq2= "<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     k x q1";
		getch();
		cout<<"\n\n\t\t\t     ("<<F1<<") x ("<<r1<<sq<<')';
		cout<<"\n\t\t\tq2= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t   ("<<k<<") x ("<<q1<<')';  getch();
		q2= (F1*r1*r1)/(k*q1);
		cout<<"\n The charge is:"<<q1<<" C";    getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : line_charge
//	DETAILS	      : IT SOLVES PROBLEMS RELATED
//			TO INFINITELY LONG CHARGED WIRE
//***************************************************************************
void line_charge( char find[] )
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"field")==0 )
	{  gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\n\t\t\t\t      ("<<lamda<<')';
		cout<<"\n\t\tE= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t     (2) x ("<<pi<<") x ("<<eps0<<") x ("<<r1<<')';  getch();

		E1= lamda/(2*pi*eps0*r1);
		cout<<"\n\n The field is:"<<E1<<" N_C(-1)"; getch();
	}
	else if( strcmp(find,"distance")==0 )
	{  cout<<"    lambda";
		cout<<"\n\t\t\tr= "<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     2"<<cpi<<ceps<<'E';
		getch();
		cout<<"\n\n\t\t\t\t    ("<<lamda<<')';
		cout<<"\n\t\tr= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t    (2) x ("<<pi<<") x ("<<eps0<<") x ("<<E1<<')'; getch();

		r1= lamda/(2*pi*eps0*E1);
		cout<<"\n\n The distance is:"<<r1<<" m";  getch();
	}
	else if( strcmp(find,"linear charge density")==0 )
	{  cout<<"    lambda= 2"<<cpi<<ceps<<"Er";
		getch();
		cout<<"\n\n    lambda= (2) x ("<<pi<<") x ("<<eps0<<") x ("<<E1<<") x ("<<r1<<')';  getch();
		lamda= 2*E1*pi*eps0*r1;
		cout<<"\n\n The linear charge density is:"<<lamda<<" C_m(-1)";  getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : plane_charge
//	DETAILS	      :	IT SOLVES PROBLEMS RELATED TO
//			CHARGED INFINITE PLANE SHEET
//***************************************************************************
void plane_charge( char find[] )
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"field")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t\t         ("<<sigma<<")";
		cout<<"\n\t\t\tE= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t    (2) x ("<<eps0<<")";   getch();

		E1= sigma/(2*eps0);
		cout<<"\n\n\t The field is:"<<E1<<" N_C(-1)"; getch();
	}
	else if( strcmp(find,"surface charge density")==0 )
	{  cout<<"    "<<cs<<" = 2 "<<ceps<<" E";
		getch();
		cout<<"\n\n\t\t    "<<cs<<" = (2) x ("<<eps0<<") x ("<<E1<<')'; getch();
		sigma= 2*E1*eps0;
		cout<<"\n\n The surface charge density is:"<<sigma<<" C_m(-2)"; getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : spherical_shell
//	DETAILS       : IT SOLVES PROBLEMS RELATED TO
//			CHARGED THIN SPHERICAL SHELL
//****************************************************************************
void spherical_shell( char find[] )
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"field")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t\t\t        "<<q1;
		cout<<"\n\t\tE= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t    (4) x ("<<pi<<") x ("<<eps0<<") x ("<<r1<<')'<<sq; getch();
		if(r1<re) E1=0;
		else
		E1= k*q1/(r1*r1);
		cout<<"\n\n The field is:"<<E1<<" N_C(-1)";   getch();
	}
	else if( strcmp(find,"distance")==0 )
	{  cout<<"  "<<(char)218<<"\t  "<<(char)191<<(char)171;
		cout<<"\n\t\t\t  "<<line<<" k x q "<<line;
		cout<<"\n\t\t\tr="<<line<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t  "<<line<<"   E   "<<line;
		cout<<"\n\t\t\t  "<<(char)192<<"\t  "<<(char)217;
		getch();
		cout<<"\n\n\t\t\t  "<<(char)218<<"\t\t      "<<(char)191<<(char)171;
		cout<<"\n\t\t\t  "<<line<<" ("<<k<<") x ("<<q1<<") "<<line;
		cout<<"\n\t\t\tr="<<line<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t  "<<line<<"\t  "<<E1<<"       "<<line;
		cout<<"\n\t\t\t  "<<(char)192<<"\t\t      "<<(char)217; getch();
		r1= sqrt( k*q1/E1 );
		cout<<"\n The distance is:"<<r1<<" m";  getch();
	}
	else if( strcmp(find,"charge")==0 )
	{  cout<<"    q= 4"<<cpi<<ceps<<"Er"<<sq;
		getch();
		cout<<"\n\n\t    q= (4) x ("<<pi<<") x ("<<eps0<<") x ("<<E1<<") x ("<<r1<<')'<<sq; getch();
		q1= E1*r1*r1/k;
		cout<<"\n\n The charge is:"<<q1<<" C"; getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : surface_area
//	DETAILS       : IT CALCULATES SURFACE AREAS
//			OF DIFFERENT 2D AND 3D OBJECTS
//***************************************************************************
void surface_area()
{  cout<<"\n Now use the formula:";
	if( strcmp(shape,"circle")==0 )
	{  cout<<"      S= "<<cpi<<'r'<<sq;
		getch();
		cout<<"\n\n\t      S= ("<<pi<<") x ("<<r1<<')'<<sq;  getch();
		S1= pi*re*re;
	}
	else if( strcmp(shape,"square")==0 )
	{  cout<<"      S= a"<<sq; getch();
		S1= re*re;
	}
	else if( strcmp(shape,"cube")==0 )
	{  cout<<"      S= 6a"<<sq; getch();
		S1= 6*re*re;
	}
	else if( strcmp(shape,"sphere")==0 )
	{  cout<<"      S= 4"<<cpi<<'r'<<sq; getch();
		S1= 4*pi*re*re;
	}
	cout<<"\n The area is:"<<S1<<" m(2)"; getch();
}
//***************************************************************************
//	FUNCTION NAME : dipole_torque
//	DETAILS       : IT SOLVES PROBLEMS RELATED TO
//			TORQUE ACTING ON A DIPOLE
//***************************************************************************
void dipole_torque( char find[] )
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"torque")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t\t   "<<ctau<<"=("<<p<<") x ("<<E1<<") x sin("<<angle<<")" ; getch();
		torque= p*E1*sin(angle);
		cout<<"\n\n The torque is:"<<torque<<" N_m";  getch();
	}
	else if( strcmp(find,"angle")==0 )
	{  cout<<"        "<<(char)218<<"     "<<(char)191;
		cout<<"\n\t\t\t        "<<line<<"  "<<ctau<<"  "<<line;
		cout<<"\n\t\t\t"<<ctheta<<" = asin"<<line<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t        "<<line<<" p E "<<line;
		cout<<"\n\t\t\t        "<<(char)192<<"     "<<(char)217;
		getch();
		cout<<"\n\n\t\t\t        "<<(char)218<<"\t\t   "<<(char)191;
		cout<<"\n\t\t\t        "<<line<<"       "<<torque<<"\t   "<<line;
		cout<<"\n\t\t\t"<<ctheta<<" = asin"<<line<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t        "<<line<<"("<<p<<") x ("<<E1<<") "<<line;
		cout<<"\n\t\t\t        "<<(char)192<<"\t\t   "<<(char)217; getch();
		angle= asin( torque/(p*E1) );
		cout<<"\n The angle is:"<<angle<<" rad"; getch();
	}
	else if( strcmp(find,"dipole moment")==0 )
	{  cout<<"      "<<ctau;
		cout<<"\n\t\t\t p= "<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    Esin"<<ctheta;
		getch();
		cout<<"\n\n\t\t\t\t   "<<torque;
		cout<<"\n\t\t\t p= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    ("<<E1<<") x sin("<<angle<<')'; getch();
		p= torque/(E1*sin(angle));
		cout<<"\n\n The dipole moment is:"<<p<<" C_m"; getch();
	}
	else if( strcmp(find,"field")==0 )
	{  cout<<"      "<<ctau;
		cout<<"\n\t\t\t E= "<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    psin"<<ctheta;
		getch();
		cout<<"\n\n\t\t\t          "<<torque;
		cout<<"\n\t\t\t E= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    ("<<p<<") x sin("<<angle<<')';  getch();
		E1= torque/(p*sin(angle));
		cout<<"\n\n The field is:"<<E1<<" N_C(-1)"; getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : gauss_law
//	DETAILS	      : IT USES GAUSS'S LAW TO CALCULATE
//			ELECTRIC FLUX, FIELD, CHARGE ETC.
//***************************************************************************
void gauss_law( char find[] )
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"flux")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t\t\t"<<q1;
		cout<<"\n\t\t\t"<<cphi<<" = "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     ("<<eps0<<')';   getch();
			flux= q1/eps0;
		cout<<"\n\n The flux is:"<<flux<<" N_m(2)_C(-1)"; getch();
	}
	else if( strcmp(find,"charge")==0 )
	{  cout<<"     q= "<<cphi<<" x "<<ceps;
		getch();
		cout<<"\n\n\t\t     q= ("<<flux<<") x ("<<eps0<<')'; getch();
		q1= flux*eps0;
		cout<<"\n\n The charge is:"<<q1<<" C";  getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : point_field
//	DETAILS		  : IT HANDLES OPREATIONS DEALING
//						 WITH FIELD PRODUCED BY A POINT CHARGE
//***************************************************************************
void point_field( char find[] )
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"field")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t\t    ("<<k<<") x ("<<q1<<')';
		cout<<"\n\t\t\tE= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t       ("<<r1<<')'<<sq; getch();
			E1=k*q1/(r1*r1);
		cout<<"\n\n The field is:"<<E1<<" N_C(-1)"; getch();
	}
	else if( strcmp(find,"charge")==0 )
	{  cout<<"        q= 4"<<cpi<<ceps<<"Er"<<sq;
		getch();
		cout<<"\n\n              q= (4) x ("<<pi<<") x ("<<eps0<<") x ("<<E1<<") x ("<<r1<<")"<<sq; getch();
		q1=E1*r1*r1/k;
		cout<<"\n\n The charge is:"<<q1<<" C";  getch();
	}
	else if( strcmp(find,"distance")==0 )
	{  cout<<"  "<<(char)218<<"      "<<(char)191<<(char)171;
		cout<<"\n\t\t\t  "<<line<<"   q  "<<line;
		cout<<"\n\t\t\tr="<<line<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t  "<<line<<" 4"<<cpi<<ceps<<"E "<<line;
		cout<<"\n\t\t\t  "<<(char)192<<"\t "<<(char)217;
		getch();
		cout<<"\n\n\t  "<<(char)218<<"\t\t\t\t\t           "<<(char)191<<(char)171;
		cout<<"\n\t  "<<line<<"\t\t\t("<<q1<<")\t\t\t   "<<line;
		cout<<"\n\tr="<<line<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t  "<<line<<" (4) x ("<<pi<<") x ("<<eps0<<") x ("<<E1<<")   "<<line;
		cout<<"\n\t  "<<(char)192<<"\t\t\t\t\t           "<<(char)217;  getch();
		r1= sqrt( k*q1/E1 );
		cout<<"\n\n The distance is:"<<r1<<" m"; getch();

	}
}
//***************************************************************************
//	FUNCTION NAME : dipole_field
//	DETAILS		  : IT HANDLES PROBLEMS RELATED
//						 TO THE FIELD PRODUCED BY A DIPOLE
//***************************************************************************
void dipole_field(char pos[])
{  cout<<"\n Now use the formula:\n\n\t\t\t";
	if( strcmp(pos,"axial")==0 )
	{  cout<<"      1      2pr ";
		cout<<"\n\t\t\t E= "<<d<<d<<d<<d<<d<<' '<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     4"<<cpi<<ceps<<"  (r"<<sq<<"- a"<<sq<<")"<<sq;
		getch();
		cout<<"\n\n\t\t\t        2 x ("<<p<<") x ("<<r1<<")";
		cout<<"\n\t\t E= ("<<k<<") x "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t{("<<r1<<')'<<sq<<"-("<<l/2<<")"<<sq<<"}"<<sq; getch();
		E1=2*k*p*r1/( pow(r1*r1-(l/2)*(l/2),2) );
	}
	else if( strcmp(pos,"equatorial")==0 )
	{  cout<<"      1        p";
		cout<<"\n\t\t\t E= "<<d<<d<<d<<d<<d<<' '<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     4"<<cpi<<ceps<<"  (r"<<sq<<"+ a"<<sq<<")^ 3/2";
		getch();
		cout<<"\n\n\t\t\t\t       ("<<p<<')';
		cout<<"\n\t\t E= ("<<k<<") x "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t       {("<<r1<<')'<<sq<<"+ ("<<l/2<<")"<<sq<<"}^ 3/2"; getch();
		E1= p*k/( pow(r1*r1+(l/2)*(l/2),3/2) );
	}
	cout<<"\n\n The field is:"<<E1<<" N_C(-1)"; getch();
}
//***************************************************************************
//	FUNCTION NAME : quantisation_of_charge
// DETAILS		  : IT USES THE FORMULA q=ne
//***************************************************************************
void quantisation_of_charge(char find[])
{  cout<<"\n Now use the given form of the formula:\n\t\t\t";
	if( strcmp(find,"number")==0 )
	{  cout<<"     q";
		cout<<"\n\t\t\t n= "<<d<<d<<d;
		cout<<"\n\t\t\t     e";
		getch();
		cout<<"\n\n\t\t\t    "<<q1;
		cout<<"\n\t\t       n= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t   "<<e; getch();
		n= q1/e;
		cout<<"\n\n The number is:"<<n;  getch();
	}
	else if( strcmp(find,"charge")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t      q = ("<<n<<") x ("<<e<<')'; getch();
		q1= n*e;
		cout<<"\n\n The charge is:"<<q1<<" C";  getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : calc_flux
//	DETAILS		  : IT CALCULATES FLUX
//***************************************************************************
void calc_flux(char find[])
{  cout<<"\n Now use the given form of the formula:\n\t\t\t";
	if( strcmp(find,"flux")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t   "<<cphi<<" = ("<<E1<<") x ("<<S1<<") x cos("<<angle<<')';  getch();
			flux= E1*S1*cos(angle);
		cout<<"\n\n The flux is:"<<flux<<" N_m(2)_C(-1)";  getch();
	}
	else if( strcmp(find,"field")==0 )
	{  cout<<"      "<<cphi;
		cout<<"\n\t\t\t E= "<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    Scos"<<ctheta;
		getch();
		cout<<"\n\n\t\t\t   ("<<flux<<')';
		cout<<"\n\t\t E= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t     ("<<S1<<") x cos("<<angle<<')'; getch();
		E1= flux/(S1*cos(angle));
		cout<<"\n\n The field is:"<<E1<<" N_C(-1)";     getch();
	}
	else if( strcmp(find,"area")==0 )
	{  cout<<"      "<<cphi;
		cout<<"\n\t\t\t S= "<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    Ecos"<<ctheta;
		getch();
		cout<<"\n\n\t\t\t   ("<<flux<<')';
		cout<<"\n\t\t S= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t     ("<<E1<<") x cos("<<angle<<')'; getch();
		S1= flux/(E1*cos(angle));
		cout<<"\n\n The area is:"<<S1<<" m(2)";         getch();
	}
	else if( strcmp(find,"angle")==0 )
	{  cout<<"\t\b"<<(char)218<<"     "<<(char)191;
		cout<<"\n\t\t\t       "<<line<<"  "<<cphi<<"  "<<line;
		cout<<"\n\t\t\t "<<ctheta<<"=acos"<<line<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t       "<<line<<" E"<<cdell<<"S "<<line;
		cout<<"\n\t\t\t\t\b"<<(char)192<<"     "<<(char)217;
		getch();
		cout<<"\n\n\t\t\t\t\b"<<(char)218<<"\t\t     "<<(char)191;
		cout<<"\n\t\t\t       "<<line<<"       ("<<flux<<")       "<<line;
		cout<<"\n\t\t\t "<<ctheta<<"=acos"<<line<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t       "<<line<<"  ("<<E1<<") x ("<<S1<<")  "<<line;
		cout<<"\n\t\t\t\t\b"<<(char)192<<"\t\t     "<<(char)217;
		angle= acos(flux/(E1*S1));
		cout<<"\n The angle is:"<<angle<<" rad"; getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : calc_dipole_moment
//	DETAILS 		  : IT CALCULATES DIPOLE MOMENT
//***************************************************************************
void calc_dipole_moment(char find[])
{   cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"dipole moment")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t       p = ("<<q1<<") x ("<<l<<")"; getch();
		p= q1*l;
		cout<<"\n\n The dipole moment is:"<<p<<" C_m";  getch();
	}
	else if( strcmp(find,"charge")==0 )
	{  cout<<"     p";
		cout<<"\n\t\t\t q= "<<d<<d<<d;
		cout<<"\n\t\t\t     l";
		getch();
		cout<<"\n\n\t\t\t     ("<<p<<")";
		cout<<"\n\t\t\t q= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     ("<<l<<")"; getch();
		q1= p/l;
		cout<<"\n The charge is:"<<q1<<" C";  getch();
	}
	else if( strcmp(find,"length")==0 )
	{  cout<<"     p";
		cout<<"\n\t\t\t l= "<<d<<d<<d;
		cout<<"\n\t\t\t     q";
		getch();
		cout<<"\n\n\t\t\t     ("<<p<<")";
		cout<<"\n\t\t\t l= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     ("<<q1<<")";	getch();
		l= p/q1;
		cout<<"\n The length is:"<<l<<" m";    getch();
	}
}
//***************************************************************************
// FUNCTION NAME : calc_linear_density
//	DETAILS		  : IT CALCULATES THE LINEAR CHARGE DENSITY
//***************************************************************************
void calc_linear_density(char find[])
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"linear charge density")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t\t          ("<<q1<<")";
		cout<<"\n\t\t\t lambda= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t  ("<<l<<")";    getch();
			lamda= q1/l;
		cout<<"\n\n The linear charge density is:"<<lamda<<" C_m(-1)"; getch();
	}
	else if( strcmp(find,"charge")==0 )
	{  cout<<"     q= lambda x l";
		getch();
		cout<<"\n\n\t\t     q= ("<<lamda<<") x ("<<l<<")"; getch();
		q1=lamda*l;
		cout<<"\n\n The charge is:"<<q1<<" C"; getch();
	}
	else if( strcmp(find,"length")==0 )
	{  cout<<"      q";
		cout<<"\n\t\t\t l= "<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    lambda";
		getch();
		cout<<"\n\n\t\t\t      ("<<q1<<")";
		cout<<"\n\t\t\t l= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t      ("<<lamda<<")";  getch();
		l= q1/lamda;
		cout<<"\n\n The length is:"<<l<<" m";  getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : calc_surface_density
//	DETAILS		  : IT CALCULATES THE SURFACE CHARGE DENSITY
//***************************************************************************
void calc_surface_density(char find[])
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"surface charge density")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t            ("<<q1<<")";
		cout<<"\n\t\t\t "<<cs<<"="<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    ("<<S1<<")";     getch();
			sigma= q1/S1;
		cout<<"\n\n The surface charge density is:"<<sigma<<" C_m(-2)"; getch();
	}
	else if( strcmp(find,"charge")==0 )
	{  cout<<"     q= "<<cs<<" x S";
		getch();
		cout<<"\n\n\t\t     q= ("<<sigma<<") x ("<<S1<<")"; getch();
		q1= sigma*S1;
		cout<<"\n\n The charge is:"<<q1<<" C";  getch();
	}
	else if( strcmp(find,"area")==0 )
	{  cout<<"     q";
		cout<<"\n\t\t\t S= "<<d<<d<<d;
		cout<<"\n\t\t\t     "<<cs;
		getch();
		cout<<"\n\n\t\t             ("<<q1<<")";
		cout<<"\n\t\t\t S= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     ("<<sigma<<")";  getch();
		S1= q1/sigma;
		cout<<"\n The area is:"<<S1<<" m(2)";  getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : calc_volume_density
//	DETAILS		  : IT CALCULATES VOLUME CHARGE DENSITY
//***************************************************************************
void calc_volume_density(char find[])
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"volume charge density")==0 )
	{   gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t\t       ("<<q1<<")";
		cout<<"\n\t\t\t rho= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t       ("<<vol<<")";   getch();
			rho= q1/vol;
		cout<<"\n\n The volume charge density is:"<<rho<<" C_m(-3)"; getch();
	}
	else if( strcmp(find,"charge")==0 )
	{  cout<<"   q= rho x V"; getch();
		cout<<"\n\n\t\t   q= ("<<rho<<") x ("<<vol<<")"; getch();
		q1= rho*vol;
		cout<<"\n\n The charge is:"<<q1<<" C";   getch();
	}
	else if( strcmp(find,"volume")==0 )
	{  cout<<"     q";
		cout<<"\n\t\t\t V= "<<d<<d<<d;
		cout<<"\n\t\t\t    rho";
		getch();
		cout<<"\n\n\t\t\t     ("<<q1<<")";
		cout<<"\n\t\t\t V= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     ("<<rho<<")";   getch();
		vol= q1/rho;
		cout<<"\n\n The volume is:"<<vol<<" m(3)";  getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : force_field
//	DETAILS		  : IT CALCULATES FORCE PRODUCED BY A FIELD
//***************************************************************************
void force_field(char find[])
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"force")==0 )
	{	 gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol();
		getch();
		cout<<"\n\n\t\t\t   F= ("<<q1<<") x ("<<E1<<")";  getch();
			F1= q1*E1;
		cout<<"\n\n The force is:"<<F1<<" N";  getch();
	}
	else if ( strcmp(find,"charge")==0 )
	{  cout<<"    F";
		cout<<"\n\t\t\t q="<<d<<d<<d;
		cout<<"\n\t\t\t    E";
		getch();
		cout<<"\n\n\t\t\t    ("<<F1<<")";
		cout<<"\n\t\t\t q="<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t    ("<<E1<<")"; getch();
		q1= F1/E1;
		cout<<"\n\n The charge is:"<<q1<<" C";   getch();
	}
	else if( strcmp(find,"field")==0 )
	{  cout<<"     F";
		cout<<"\n\t\t\t E= "<<d<<d<<d;
		cout<<"\n\t\t\t     q";
		getch();
		cout<<"\n\n\t\t\t     ("<<F1<<")";
		cout<<"\n\t\t\t E= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t     ("<<q1<<")";  getch();
		E1= F1/q1;
		cout<<"\n\n The field is:"<<E1<<" N_C(-1)";   getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : calc_volume
//	DETAILS		  : IT CALCULATES VOLUME OF DIFFERENT SHAPES
//***************************************************************************
void calc_volume()
{	if( strcmp(shape,"sphere")==0 )
	{  cout<<"  V = 4/3 "<<cpi<<"r^3"; getch();
		vol= 4/3*pi*re*re*re;
	}
	else if( strcmp(shape,"cube")==0 )
	{  cout<<"    V= r^3"; getch();
		vol= re*re*re;
	}
	cout<<"\n The volume is:"<<vol<<" m(3)"; getch();
}
//***************************************************************************
//	FUNCTION NAME : cell_density
//	DETAILS		  : IT SOLVES PROBLEMS RELATED TO THE FORMULA
//						 FOR CALCULATION OF CELL DENSITY
//***************************************************************************
void cell_density(char find[])
{  cout<<"\n Now use the given form of the formula:\n\n\t\t\t";
	if( strcmp(find,"density")==0 )
	{  gotoxy(1,wherey()-2);
		clreol();
		cout<<"\n Use the formula in the same form";clreol(); getch();
		density=z*M/(a*a*a*N_a);
		cout<<"\n\n The density is:"<<density; getch();
	}
	else if( strcmp(find,"molar mass")==0 )
	{  cout<<"     d x a^3 x N_a";
		cout<<"\n\t\t\t  M= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t  z";
		getch();
		cout<<"\n\n\t\t\t     ("<<density<<") x ("<<a*a*a<<") x ("<<N_a<<')';
		cout<<"\n\t\t\t  M= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t\t  ("<<z<<')';  getch();
		M=density*a*a*a*N_a/z;
		cout<<"\n\n The molar mass is:"<<M<<" g_mol(-1)"; getch();
	}
	else if( strcmp(find,"nature")==0 )
	{  cout<<"     d x a^3 x N_a";
		cout<<"\n\t\t\t  z= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t  M";
		getch();
		cout<<"\n\n\t\t\t     ("<<density<<") x ("<<a*a*a<<") x ("<<N_a<<')';
		cout<<"\n\t\t\t  z= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t\t  ("<<M<<')';  getch();
		z=density*a*a*a*N_a/M;
		if(z==1)
		cout<<"\n\n The unit cell is in simple packing.";
		else if(z==2)
		cout<<"\n\n The unit cell is in body centred cubic packing.";
		else if(z==4)
		cout<<"\n\n The unit cell is in face centred cubic packing (/ccp/hcp).";
		getch();
	}
	else if( strcmp(find,"edge")==0 )
	{	cout<<"     "<<(char)218<<"\t\t\b"<<(char)191<<"1/3";
		cout<<"\n\t\t\t     "<<line<<"  z x M  "<<line;
		cout<<"\n\t\t\t  a= "<<line<<d<<d<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t     "<<line<<"  d x N_a "<<line;
		cout<<"\n\t\t\t     "<<(char)192<<"         "<<(char)217;
		getch();
		cout<<"\n\t\t\t     "<<(char)218<<"\t\t\t\t\b\b"<<(char)191<<"1/3";
		cout<<"\n\t\t\t     "<<line<<"  ("<<z<<" x ("<<M<<')'<<"\t\b\b"<<line;
		cout<<"\n\t\t\t  a= "<<line<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<line;
		cout<<"\n\t\t\t     "<<line<<"  ("<<density<<") x ("<<N_a<<") "<<line;
		cout<<"\n\t\t\t     "<<(char)192<<"      \t\t      "<<(char)217;
		getch();
		a= pow( z*M/(density*N_a),1.0/3.0);
		cout<<"\n\n The edge length is:"<<a<<" m";
		getch();
	}
}
//***************************************************************************
//	FUNCTION NAME : radius_edge
//	DETAILS 		  : IT CALCULATES RADIUS FROM CELL EDGE
//***************************************************************************
void radius_edge(char find[])
{  cout<<"\n\n Now use the formula:\n\t\t\t";
	if( strcmp(find,"radius")==0 )
	{	if( strcmp(pack,"simple")==0 )
		{  cout<<"     a";
			cout<<"\n\t\t\t r= "<<d<<d<<d;
			cout<<"\n\t\t\t     2";
			getch();
			cout<<"\n\n\t\t\t    ("<<a<<')';
			cout<<"\n\t\t\t r= "<<d<<d<<d<<d<<d<<d<<d;
			cout<<"\n\t\t\t       2";
			re=a/2;
			cout<<"\n\n The radius is:"<<re<<" m";
			getch();
		}
		else if( strcmp(pack,"bcc")==0 )
		{  cout<<"     "<<sqr<<"3";
			cout<<"\n\t\t\t  r= "<<d<<d<<d<<" a";
			cout<<"\n\t\t\t      2";
			getch();
			cout<<"\n\t\t\t     "<<sqr<<"3";
			cout<<"\n\t\t\t  r= "<<d<<d<<d<<" ("<<a<<')';
			cout<<"\n\t\t\t      2";
			getch();
			re= sqrt(3)*a/4;
			cout<<"\n\n The radius is:"<<re<<" m";
			getch();
		}
		else if( strcmp(pack,"fcc")==0 )
		{  cout<<"      a";
			cout<<"\n\t\t\t  r= "<<d<<d<<d;
			cout<<"\n\t\t\t     2"<<sqr<<"2";
			getch();
			cout<<"\n\n\t\t\t     ("<<a<<')';
			cout<<"\n\t\t\t  r= "<<d<<d<<d<<d<<d<<d<<d;
			cout<<"\n\t\t\t       2"<<sqr<<"2";
			getch();
			re= a/( sqrt(2)*2 );
			cout<<"\n\n The radius is:"<<re<<" m";
			getch();
		}
	}
	else if( strcmp(find,"edge")==0 )
	{	if( strcmp(pack,"simple")==0 )
		{  cout<<"      a= 2 x r";
			getch();
			cout<<"\n\n\t\t\t   a= 2 x ("<<re<<')';
			getch();
			a=2*re;
			cout<<"\n\n The edge length is:"<<a<<" m";
			getch();
		}
		else if( strcmp(pack,"bcc")==0 )
		{  cout<<"       4";
			cout<<"\n\t\t\t   a= "<<d<<d<<d<<" r";
			cout<<"\n\t\t\t      "<<sqr<<'3';
			getch();
			cout<<"\n\n\t\t\t       4";
			cout<<"\n\t\t\t   a= "<<d<<d<<d<<" ("<<re<<')';
			cout<<"\n\t\t\t      "<<sqr<<'3';
			getch();
			a= 4*re/sqrt(3);
			cout<<"\n\n The edge length is:"<<a<<" m";
			getch();
		}
		else if( strcmp(pack,"fcc")==0 )
		{  cout<<"        a= 2"<<sqr<<"2r";
			getch();
			cout<<"\n\n\t\t\t   a= 2 x"<<sqr<<"2 x ("<<re<<')';
			getch();
			a= sqrt(2)*2*re;
			cout<<"\n\n The edge length is:"<<a<<" m";
			getch();
		}
	}
}

//***************************************************************************
//	FUNCTION NAME : cell_density
// DETAILS		  :
//***************************************************************************
/*void cell_density(char find[])
{	if( strcmp(find,"density")==0 )
	{  density= z*M/(a*a*a*N_a);
	}
	else if( strcmp(find,"molar mass")==0 )
	{	M= density*N_a*a*a*a/z;
	}
	else if( strcmp(find,"edge")==0 )
	{	a= pow( z*M/(density*N_a),1.0/3.0);
	}
	else if( strcmp(find,"nature")==0 )
	{	z= density*N_a*a*a*a/M;
	}
}
//***************************************************************************
//	FUNCTION NAME : radius_edge
//	DETAILS		  :
//***************************************************************************
void radius_edge(char find[])
{	if( strcmp(find,"radius")==0 )
	{	if( strcmp(pack,"simple")==0 )
		{	re=a/2;
		}
		else if( strcmp(pack,"bcc")==0 )
		{	re= sqrt(3)/4.0*a;
		}
		else if( strcmp(pack,"fcc")==0 )
		{	re= a/( 2*sqrt(2) );
		}
	}
	else if( strcmp(find,"edge")==0 )
	{	if( strcmp(pack,"simple")==0 )
		{	a= 2*re;
		}
		else if( strcmp(pack,"bcc")==0 )
		{	a= 4.0/sqrt(3)*re;
		}
		else if( strcmp(pack,"fcc")==0 )
		{	a= re*2*sqrt(2);
		}
	}
} */
//###########################################################################
//###########################################################################
//***************************************************************************
//	FUNCTION NAME : check_dipole_length
//	DETAILS		  : IT CHECKS WHETHER DIPOLE LENGTH IS GIVEN
//***************************************************************************\
void check_dipole_length()
{	for(int i=0;i<5;++i)
	{	if( to_find[i]==3 )
		{	if( qsign.sign[5]==1 )
			{	l=r1; qsign.sign[13]++; break; }
		}
	}
}
//***************************************************************************
//	FUNCTION NAME : data_initalizer
//	DETAILS       : IT TAKES THE NUMERICAL DATA
//			ALONG WITH UNITS AND INITIALIZES
//			THE DESIRED LOGICAL VARIABLES
//***************************************************************************
void data_initializer_phy( char variable[])
{	if( strcmp(variable,"force")==0 )
		F[qsign.sign[0]++] = data.retmag();
	else if( strcmp(variable,"charge")==0 )
		q[qsign.sign[1]++] = data.retmag();
	else if( strcmp(variable,"area")==0 )
		S[qsign.sign[9]++] = data.retmag();
	else if ( strcmp(variable,"distance")==0 )
		r[qsign.sign[5]++]= data.retmag();
	else if( strcmp(variable,"field")==0 )
		E[qsign.sign[2]++]= data.retmag();
	else if( strcmp(variable,"dipole moment")==0 )
	{   p= data.retmag(); qsign.sign[3]++; }
	else if( strcmp(variable,"flux")==0 )
	{	flux= data.retmag(); qsign.sign[4]++; }
	else if( strcmp(variable,"linear charge density")==0 )
	{	lamda= data.retmag(); qsign.sign[6]++; }
	else if( strcmp(variable,"surface charge density")==0 )
	{	sigma= data.retmag(); qsign.sign[7]++; }
	else if( strcmp(variable,"volume charge density")==0 )
	{	rho= data.retmag(); qsign.sign[8]++; }
	else if(strcmp(variable,"torque")==0 )
	{	torque= data.retmag(); qsign.sign[10]++; }
	else if( strcmp(variable,"angle")==0 )
	{	angle= data.retmag(); qsign.sign[11]++; }
	else if( strcmp(variable,"number")==0 )
	{	n= data.retmag();	qsign.sign[12]++; }
	else if( strcmp(variable,"radius")==0 )
	{	re= data.retmag(); qsign.sign[14]++; }
	else if( strcmp(variable,"diameter")==0 )
	{	re =data.retmag()/2; qsign.sign[14]++; }
	else if( strcmp(variable,"edge")==0 )
	{	re= data.retmag(); qsign.sign[14]++; }
	else if( strcmp(variable,"volume")==0 )
	{	vol= data.retmag(); qsign.sign[17]++; }
}
//***************************************************************************
//	FUNCTION NAME : data_extracter
//	DETAILS	      : IT EXTRACTS NUMERICAL DATA
//			ALONG WITH UNITS AND CALL data_initializer_phy
//***************************************************************************
void data_extracter_phy()
{
	for(int i=0,j,flag;i<tw;++i)
	{
	if( strcmp(input[i],"axis")==0||strcmp(input[i],"axial")==0 )
		{	strcpy(position,"axial"); qsign.sign[15]++; }
		if( strcmp(input[i],"equatorial")==0 )
		{	strcpy(position,"equatorial"); qsign.sign[15]++; }

		if( isdigit(input[i][0]) || isdigit(input[i][1]) )
		{
			data.getdata( atof(input[i]), input[i+1] );
			convert_to_SI(data);
			if( strcmp(data.retunit(),"N")==0 )
				data_initializer_phy("force");
			else if( strcmp(data.retunit(),"C")==0 )
				data_initializer_phy("charge");
			else if( strcmp(data.retunit(),"m(2)")==0 )
				data_initializer_phy("area");
			else if( strcmp(data.retunit(),"m(3)")==0 )
				data_initializer_phy("volume");
			else if( strcmp(data.retunit(),"V_m(-1)")==0  )
				data_initializer_phy("field");
			else if( strcmp(data.retunit(),"N_C(-1)")==0 )
				data_initializer_phy("field");
			else if( strcmp(data.retunit(),"C_m")==0 )
				data_initializer_phy("dipole moment");
			else if( strcmp(data.retunit(),"V_m")==0 )
				data_initializer_phy("flux");
			else if( strcmp(data.retunit(),"N_m(2)_C(-1)")==0 )
				data_initializer_phy("flux");
			else if( strcmp(data.retunit(),"C_m(-1)")==0 )
				data_initializer_phy("linear charge density");
			else if( strcmp(data.retunit(),"C_m(-2)")==0 )
				data_initializer_phy("surface charge density");
			else if( strcmp(data.retunit(),"C_m(-3)")==0 )
				data_initializer_phy("volume charge density");
			else if( strcmp(data.retunit(),"N_m")==0 )
				data_initializer_phy("torque");
			else if( strcmp(data.retunit(),"rad")==0 )
				data_initializer_phy("angle");
			else if( strcmp(data.retunit(),"m")==0 )
			{	for(j=-6,flag=0;j<=6;++j)
				{	if( strcmp(input[i+j],"radius")==0 )
					{	data_initializer_phy("radius"); flag=1;  break; }
					else if( strcmp(input[i+j],"diameter")==0 || strcmp(input[i+j],"diametre")==0 )
					{	data_initializer_phy("diameter"); flag=1; break; }
					else if( strcmp(input[i+j],"side")==0 || strcmp(input[i+j],"edge")==0 )
					{	data_initializer_phy("edge"); flag=1; break; }
				}
				if(!flag)
					data_initializer_phy("distance");
			}
			else
			data_initializer_phy("number");
		}
	}
}
//***************************************************************************
//	FUNCTION NAME : data_extracter_chem
//	DETAILS		  :
//***************************************************************************
void data_extracter_chem()
{	for(int j,i=0;i<tw;++i)
	{	if( isdigit(input[i][0]) || isdigit(input[i][1]) )
		{
			data.getdata( atof(input[i]), input[i+1] );
			convert_to_SI(data);
			if( strcmp(data.retunit(),"g_mol")==0 )
			  data_initializer_chem("molar mass");
			else if( strcmp(data.retunit(),"kg_m(-3)")==0 )
				data_initializer_chem("density");
			else if( strcmp(data.retunit(),"m")==0 )
			{	for(j=-6;j<=6;j++)
				{  if( strcmp(input[i+j],"radius")==0 )
					{	data_initializer_chem("radius");
						break;
					}
					else if( strcmp(input[i+j],"diameter")==0 || strcmp(input[i+j],"diametre")==0 )
					{	data_initializer_chem("diameter");
						break;
					}
					else if( strcmp(input[i+j],"side")==0 || strcmp(input[i+j],"edge")==0 )
					{	data_initializer_chem("edge");
						break;
					}
				}
			}
			else	data_initializer_chem("number");
		}
	}
}
//***************************************************************************
//	FUNCTION NAME : data_initializer_chem
//	DETAILS		  :
//***************************************************************************
void data_initializer_chem(char variable[])
{  if( strcmp(variable,"molar mass")==0 )
	{  M= data.retmag(); qsign.sign[0]++; }
	else if( strcmp(variable,"density")==0 )
	{	density= data.retmag(); qsign.sign[4]++; }
	else if( strcmp(variable,"radius")==0 )
	{  re= data.retmag(); qsign.sign[2]++; }
	else if( strcmp(variable,"diameter")==0 )
	{	re= data.retmag()/2; qsign.sign[2]++; }
	else if( strcmp(variable,"edge")==0 )
	{	a= data.retmag(); qsign.sign[1]++; }
	else if( strcmp(variable,"number")==0 )
	{	z= data.retmag(); qsign.sign[3]++; }
}
//***************************************************************************
//	FUNCTION NAME : shape_identifier
//	DETAILS		  : IT IDENTIFIES THE GIVEN 2D OR 3D
//						 OBJECT IN THE PROBLEM
//***************************************************************************
void shape_identifier()
{	for(int i=0;i<tw;++i)
	{	if( strcmp(input[i],"cube")==0 || strcmp(input[i],"cubic")==0 )
		{	strcpy(shape,"cube");qsign.sign[16]++; break; }
		else if( strcmp(input[i],"sphere")==0 || strcmp(input[i],"spherical")==0 )
		{	strcpy(shape,"sphere");qsign.sign[16]++; break; }
		else if( strcmp(input[i],"square")==0 )
		{	strcpy(shape,"square");qsign.sign[16]++; break; }
		else if( strcmp(input[i],"circle")==0 || strcmp(input[i],"circular")==0 )
		{	strcpy(shape,"circle");qsign.sign[16]++; break; }
	}
}
//***************************************************************************
//	FUNCTION NAME : pack_identifier
//	DETAILS		  :
//***************************************************************************
void pack_identifier()
{	for(int i=0;i<tw;++i)
	{	if( !strcmp(input[i],"fcc") || !strcmp(input[i],"hcp") || !strcmp(input[i],"ccp") )
		{	strcpy(pack,"fcc"); z=4; }
		else if( !strcmp(input[i],"bcc") )
		{	strcpy(pack,"bcc"); z=2; }
		else if( strcmp(input[i],"cubic")==0 )
		{	if(!strcmp(input[i-1],"face-centred") || !strcmp(input[i],"hexagonal") )
			{	strcpy(pack,"fcc"); z=4; }
			else if( strcmp(input[i+1],"close")==0 )
			{	strcpy(pack,"fcc"); z=4; }
			else if(strcmp(input[i-1],"body-centred")==0 )
			{	strcpy(pack,"bcc"); z=2; }
			else if( !strcmp(input[i-1],"simple") )
			{	strcpy(pack,"simple"); z=1; }
		}
	}
}
//***************************************************************************
//	FUNCTION NAME : aim
//	DETAILS 		  : IT SEARCHES FOR WHAT IS TO BE CALCULATED
//***************************************************************************
void aim_phy()
{  int flag;
	for(int i=0,j;i<tw;++i)
	{	if( aim_phy(input[i]) )
		{	for(j=1;j<11;j++)
			{	aim_phy(input[i+j],flag,i+j);
				if(flag)
				{	if( strcmp(input[i+j+1],"and")==0 )
						aim_phy(input[i+j+2],flag,i+j+2);
				}
			}
		}
	}
}

int aim_phy(char str[])
{	ifstream fin("B_QWORDS.dat",ios::in|ios::binary);
	words in;
	if(!fin)
	{	cout<<"\n B_QWORDS.dat not found";
		getch();
	}
	fin.seekg(0);
	while(fin)
	{  fin.read( (char*)&in,sizeof(in) );
		if( strcmp(str,in.keyword)==0 )
		{ fin.close(); return 1; }
	}
	fin.close();
	return 0;
}

void aim_phy(char str[],int& f,int i)
{	ifstream fin("B_PQWORD.dat",ios::in|ios::binary);
	words in;
	if(!fin)
	{	cout<<"\n B_PQWORD.dat not found";
		getch();
	}
	fin.seekg(0);
	while(fin)
	{  fin.read( (char*)&in,sizeof(in) );
		if( strcmp(str,in.keyword)==0 )
		{	to_find[itf++]=case_selecter_phy(str);
			qsign.sign[case_selecter_phy(str)]++;
			f=1;
		}
		else if( strcmp(str,"charge")==0 )
		{	if( !strcmp(input[i-1],"linear") && !strcmp(input[i+1],"density") )
			{	to_find[itf++]=6;
				qsign.sign[6]++; f=1;break;
			}
			else if(!strcmp(input[i-1],"surface") && !strcmp(input[i+1],"density") )
			{	to_find[itf++]=7;
				qsign.sign[7]++;f=1; break;
			}
			else if(!strcmp(input[i-1],"volume") && !strcmp(input[i+1],"density") )
			{	to_find[itf++]=8;
				qsign.sign[8]++; f=1;break;
			}
			else
			{	to_find[itf++]=1;
				qsign.sign[1]++; f=1;break;
			}

		}
		else if( !strcmp(str,"dipole") && !strcmp(input[i+1],"moment") )
		{	to_find[itf++]=3;
			qsign.sign[3]++;
			f=1;   break;
		}
		else if( !strcmp(str,"volume") && strcmp(input[i+1],"charge") )
		{	to_find[itf++]=17;
			qsign.sign[17]++;
			f=1; break;
		}
	}
	f=0;
}

void aim_chem()
{  int flag;
	for(int i=0,j;i<tw;++i)
	{	if( aim_chem( input[i] ) )
		{	for(j=1;j<11;j++)
			{	aim_chem(input[i+j],flag,i+j);
				if(flag)
				{	if( strcmp(input[i+j+1],"and")==0 )
						aim_chem(input[i+j+2],flag,i+j+2);
				}
			}
		}
	}
}

int aim_chem(char str[])
{	ifstream fin("B_QWORDS.dat",ios::in|ios::binary);
	words in;
	if(!fin)
	{	cout<<"\n B_QWORDS.dat not found";
		getch();
	}
	fin.seekg(0);
	while(fin)
	{  fin.read( (char*)&in,sizeof(in) );
		if( strcmp(str,in.keyword)==0 )
		{ fin.close(); return 1; }
	}
	fin.close();
	return 0;
}

void aim_chem(char str[],int& f,int i)
{	ifstream fin("B_CQWORD.dat",ios::in|ios::binary);
	words in;
	if(!fin)
	{	cout<<"\n B_CQWORD.dat not found";
		getch();
	}
	fin.seekg(0);
	while(fin)
	{  fin.read( (char*)&in,sizeof(in) );
		if( strcmp(str,in.keyword)==0 )
		{	to_find[itf++]=case_selecter_chem(str);
			qsign.sign[case_selecter_chem(str)]++;
			f=1;
		}
		else if( strcmp(str,"mass")==0 )
		{	if( !strcmp(input[i-1],"molar") || !strcmp(input[i-1],"molecular") || !strcmp(input[i-1],"atomic") )
			{	to_find[itf++]=6;
				qsign.sign[6]++; f=1;break;
			}
		}
	}
	f=0;
}

//***************************************************************************
//	FUNCTION NAME : remove_punctuation
//	DETAILS		  : IT REMOVES PUNCTUATION MARKS
//						 ATTACHED WITH THE WORDS
//***************************************************************************
void remove_punctuation(char A[])
{	int i=strlen(A)-1;
	if(A[i]==','||A[i]=='?'||A[i]=='!'||A[i]=='.'||A[i]==':'||A[i]==';'||A[i]==')'||A[i]=='"'||A[i]=='\'')
		A[i]='\0';
	if(A[0]=='('||A[0]=='\''||A[0]=='"')
		for(i=0;A[i]!='\0';++i)
			A[i]=A[i+1];
}
//***************************************************************************
//	FUCTION NAME : function_selecter
//	DETAILS 		 :	IT SELECTS THE REQUIRED FORMULA FUNCTION
//***************************************************************************
void function_selecter_phy()
{  int flag=1;
	signature in;
	ifstream fin("B_QSIGN.dat",ios::in|ios::binary);
	fin.seekg(0);
	while(fin)
	{	fin.read( (char*)&in,sizeof(in) );
		if( sign_matcher_full_phy(in) )
		{	function_caller_phy(in.n);
			flag=0;
			break;
		}

	}
	if(flag)
	{	fin.seekg(0);
		while(fin)
		{	fin.read( (char*)&in,sizeof(in) );
			if( signmp(in) )
			{	function_caller_phy(in.n);
				break;
			}
		}
	}
	fin.close();
}

void function_selecter_chem()
{  int flag=1;
	signature in;
	ifstream fin("B_CQSIGN.dat",ios::in|ios::binary);
	fin.seekg(0);
	while(fin)
	{	fin.read( (char*)&in,sizeof(in) );
		if( sign_matcher_full_chem(in) )
		{	function_caller_chem(in.n);
			flag=0;
			break;
		}

	}
	if(flag)
	{	fin.seekg(0);
		while(fin)
		{	fin.read( (char*)&in,sizeof(in) );
			if( sign_matcher_inc_chem(in) )
			{	function_caller_chem(in.n);
				break;
			}
		}
	}
	fin.close();
}
//***************************************************************************
//	FUNCTION NAME : sign_matcher
//	DETAILS		  : IT MATCHES THE SIGNATURE OF THE QUESTION
//						 WITH THAT STORED IN THE FILE
//***************************************************************************
int sign_matcher_full_phy(signature in)
{  int count=0;
	{	for(int i=0;i<18;++i)
		{	if( in.sign[i]!=0 && qsign.sign[i]==in.sign[i] )
					count++;
		}
	}
	if(count==in.total)
	return 1;
	return 0;
}

int signmp(signature in)
{  int count=0;
	if( in.sign[to_find[itf-1]] )
	{	for(int i=0;i<18;++i)
		{	if( in.sign[i]!=0  && qsign.sign[i]==in.sign[i] )
				count++;
		}
	}
	if(count==in.total-1)
	{	for(int i=0;i<18;++i)
		{  if( in.sign[i]!=0  && !qsign.sign[i]==in.sign[i] )
			{	to_find[itf++]=i;
				qsign.sign[i]++;
				function_selecter_phy();
				return 1;
			}
		}
	 }
	return 0;
}

int sign_matcher_full_chem(signature in)
{  int count=0;
	{	for(int i=0;i<5;++i)
		{	if( in.sign[i]!=0 && qsign.sign[i]==in.sign[i] )
					count++;
		}
	}
	if(count==in.total)
	return 1;
	return 0;
}

int sign_matcher_inc_chem(signature in)
{  int count=0;
	if( in.sign[to_find[itf-1]] )
	{	for(int i=0;i<5;++i)
		{	if( in.sign[i]!=0  && qsign.sign[i]==in.sign[i] )
				count++;
		}
	}
	if(count==in.total-1)
	{	for(int i=0;i<5;++i)
		{  if( in.sign[i]!=0  && !qsign.sign[i]==in.sign[i] )
			{	to_find[itf++]=i;
				qsign.sign[i]++;
				function_selecter_chem();
				return 1;
			}
		}
	 }
	return 0;
}

//***************************************************************************
//***************************************************************************
//	FUNCTION NAME : function_caller
//	DETAILS		  : IT CALLS THE REQUIRED FORMULA FUNCTION
//***************************************************************************
void function_caller_phy(int num)
{  case_selecter_phy( to_find[itf-1] );
	cout<<"\n\n Now use the formula:";
	switch(num)
	{
		case 0 : cout<<"      k x q1 x q2";
					cout<<"\n\t\t\tF= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
					cout<<"\n\t\t\t        r"<<sq;  getch();
					coloumbs_law(cstring);
					break;

		case 1 : cout<<"       lambda";
					cout<<"\n\t\t\tE= "<<d<<d<<d<<d<<d<<d<<d<<d;
					cout<<"\n\t\t\t     2"<<cpi<<ceps<<'r';  getch();
					line_charge(cstring);
					break;

		case 2 : cout<<"       "<<cs;
					cout<<"\n\t\t\tE= "<<d<<d<<d<<d;
					cout<<"\n\t\t\t    2"<<ceps;   getch();
					plane_charge(cstring);
					break;

		case 3 : cout<<"        q";
					cout<<"\n\t\t\tE= "<<d<<d<<d<<d<<d<<d;
					cout<<"\n\t\t\t    4"<<cpi<<ceps<<'r'<<sq; getch();
					spherical_shell(cstring);
					break;

		case 4 : gotoxy(1,wherey()); clreol();
					surface_area();
					break;

		case 5 : cout<<"\t   "<<carr<<"  "<<carr<<"   "<<carr;
					cout<<"\n\t\t\t   "<<ctau<<"= p x E = pEsin"<<ctheta; getch();
					dipole_torque(position);
					break;

		case 6 : gotoxy(1,wherey()); clreol();
					dipole_field(cstring);
					break;

		case 7 : cout<<"        q";
					cout<<"\n\t\t\t"<<cphi<<" = "<<d<<d<<d;
					cout<<"\n\t\t\t     "<<ceps;   getch();
					gauss_law(cstring);
					break;

		case 8 : cout<<"        q";
					cout<<"\n\t\t\tE= "<<d<<d<<d<<d<<d<<d;
					cout<<"\n\t\t\t    4"<<cpi<<ceps<<'r'<<sq; getch();
					point_field(cstring);
					break;

		case 9 : cout<<"\t     "<<"  "<<carr<<"   "<<carr;
					cout<<"\n\t\t\t   "<<cphi<<" = E x "<<cdell<<"S = E"<<cdell<<"Scos"<<ctheta;
					getch();
					calc_flux(cstring);
					break;

		case 10 :cout<<"      p = q x l"; getch();
					calc_dipole_moment(cstring);
					break;

		case 11 :cout<<"      q = n x e"; getch();
					quantisation_of_charge(cstring);
					break;

		case 12 :cout<<"             q";
					cout<<"\n\t\t\t lambda= "<<d<<d<<d;
					cout<<"\n\t\t\t\t  l";    getch();
					calc_linear_density(cstring);
					break;

		case 13 :cout<<"       q";
					cout<<"\n\t\t\t "<<cs<<"="<<d<<d<<d;
					cout<<"\n\t\t\t    S";     getch();
					calc_surface_density(cstring);
					break;

		case 14 :cout<<"          q";
					cout<<"\n\t\t\t rho= "<<d<<d<<d;
					cout<<"\n\t\t\t       V";   getch();
					calc_volume_density(cstring);
					break;

		case 15 :cout<<"\t   "<<carr<<"      "<<carr;
					cout<<"\n\t\t\t   F= q x E";  getch();
					force_field(cstring);
					break;

		case 16 :calc_volume();
					break;

		default : ;
	}
	itf--;
}

void function_caller_chem(int num)
{  case_selecter_chem( to_find[itf-1] );
	cout<<"\n\n Now use the formula:";
	switch(num)
	{
		case 0 :	cout<<"          z x M ";
					cout<<"\n\t\t\t  d= "<<d<<d<<d<<d<<d<<d<<d<<d<<d;
					cout<<"\n\t\t\t      a^3 x Na";
					getch();
					cell_density(cstring);
					break;

		case 1 :
					radius_edge(cstring);
					break;

		default : ;
	}
	itf--;
}


//***************************************************************************
//	FUNCTION NAME : summary
//	DETAILS		  : IT DISPLAYS THE SUMMARY OF THE UNIT
//***************************************************************************
void summary()
{	int n=1,key; char ch;
	ifstream fin("SUM1.txt");
	clrscr();
	while(!fin.eof())
	{	fin.get(ch);
		cout<<ch;
	}
	fin.close();
	A:
	while(bioskey(1)==0);
	key=bioskey(0);
	if(key==19200)
	{	if(n!=1) n--;
		else goto A;
	}
	else if(key==19712)
	{	if(n!=9) n++;
		else goto A;
	}
	else if(key==283) pcmenu1();
	else goto A;

	clrscr();
	if(n==1){   fin.open("SUM1.txt");
					fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==2){ fin.open("SUM2.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	 }
	else if(n==3){ fin.open("SUM3.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==4){ fin.open("SUM4.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==5){ fin.open("SUM5.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==6){ fin.open("SUM6.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==7){ fin.open("SUM7.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==8){ fin.open("SUM8.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==9){ fin.open("SUM9.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();

	}
	goto A;
}
//***************************************************************************
//	FUNCTION NAME : help
//	DETAILS		  : IT DISPLAYS THE USER MANUAL
//***************************************************************************
void help()
{	int n=1,key; char ch;
	ifstream fin("HELP1.txt");
	clrscr();
	while(!fin.eof())
	{	fin.get(ch);
		cout<<ch;
	}
	fin.close();
	A:
	while(bioskey(1)==0);
	key=bioskey(0);
	if(key==19200)
	{	if(n!=1) n--;
		else goto A;
	}
	else if(key==19712)
	{	if(n!=2) n++;
		else goto A;
	}
	else if(key==283 || key==7181) pcmenu3();
	else goto A;

	clrscr();
	if(n==1){   fin.open("HELP1.txt");
					fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==2){ fin.open("HELP2.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	 }
	 goto A;

}
//***************************************************************************
//	FUNCTION NAME : getwords
//	DETAILS	      : IT BREAKS THE INPUT STATEMENT
//						  DISCRETE WORDS
//***************************************************************************
void getwords()
{       int i,j,k;
	for(i=0,j=0,k=0;tmp[i]!='\0';++i)
	{	if(tmp[i]!=' ')
			input[k][j++]=tmp[i];
		else
		{	k++;
			j=0;
		}
	}
	tw=k+1;
	for(i=0;i<tw;++i)
	remove_punctuation(input[i]);
}
//**************************************************************************
// FUNCTION NAME : scr
//	DETAILS		  : IT PRINTS "HELLO"
//***************************************************************************
void scr()
{
clrscr();
cout<<"\n\n\n\n\n\n\n";
cout<<"\n\t\t |||    |||   ||||||\t||\t||     |||\n\t\t ||||||||||   ||||  \t||  \t||   ||   ||";

cout<<"\n\t\t |||    |||   ||||||\t||||\t||||   |||  .........";

}
//***************************************************************************
//	FUNCTION NAME : jnv
//	DETAILS		  : IT PRINTS "JNV"
//***************************************************************************
void jnv()
{
cout<<"\n \n ";
cout<<"\n\t\t\t\t||||||\t||  ||\t||  ||\n\t\t\t\t  ||  \t||\\ ||\t |  |\n\t\t\t\t||||  \t|| \\||\t  ||";
}
//***************************************************************************
//	FUNCTION NAME : name
//	DETAILS		  : IT PRINTS THE PROGRAM'S NAME
//***************************************************************************
void name()
{clrscr();
textcolor(LIGHTRED+LIGHTBLUE);
cout<<"\n\n";
cprintf("                    * * * * * *  * *     * *  * * * * * *");cout<<endl;
cprintf("                    * * * * * *  * *     * *  * * * * * *");cout<<endl;
cprintf("                        * *      * *     * *      * *");cout<<endl;
cprintf("                        * *      * * * * * *      * *");cout<<endl;
cprintf("                        * *      * * * * * *      * *");cout<<endl;
cout<<"\n\n\n\n\n";
cprintf("                    ( T H E  U L T I M A T E  T U T O R )");
cout<<"\n\n\n\n";
cprintf("            D E V E L O P E D  B Y  :-  S O U R A B H  J A I S W A L");
}

//****************************************************************************
//	FUNCTION NAME : exitdesign
//	DETAILS		  : IT PRINTS THE EXIT TIME DESIGN
//****************************************************************************
void exitdesign()
{  clrscr();
	int i=1;
	while(bioskey(1)==0)
	{ delay(200);
	textcolor(i++);
	cout<<"\n";
	cprintf("     *-*    *      *-*-*    *    *-*-*  *-*-*  *--*        *****    *");
	cout<<endl;
	delay(100);
	cprintf("     *   *  *      *       / \\   *      *      *   \\         *     * *");
	cout<<endl;
	delay(100);
	cprintf("     *-*    *      *-*    *-*-*  *-*-*  *-*    *    *        *    *   *");
	cout<<endl;delay(100);
	cprintf("     *      *      *      *   *      *  *      *   /         *     * *");
	cout<<endl;delay(100);
	cprintf("     *      *-*-*  *-*-*  *   *  *-*-*  *-*-*  *--*          *      *");
	cout<<"\n\n\n"; delay(100);
	cprintf("         *   *  *-*-*  *      *-*         *   *    *    *   *"); cout<<endl;delay(100);
	cprintf("         *   *  *      *      *   *        * *    * *   *   *"); cout<<endl;delay(100);
	cprintf("         *-*-*  *-*    *      *-*           *    *   *  *   *"); cout<<endl;delay(100);
	cprintf("         *   *  *      *      *             *     * *   *   *"); cout<<endl;delay(100);
	cprintf("         *   *  *-*-*  *-*-*  *             *      *    *-*-*");
	cout<<"\n\n\n\n";
	cprintf("   *-*-*  *-*-*  *-*-*    *   *    *    *   *     *-*-*    *      *    *   *"); cout<<endl;delay(100);
	cprintf("   *      *      *         * *    * *   *   *     *       * *    * *   **  *"); cout<<endl;delay(100);
	cprintf("   *-*-*  *-*    *-*        *    *   *  *   *     *-*-*  *   *  *   *  * * *"); cout<<endl;delay(100);
	cprintf("       *  *      *          *     * *   *   *         *   * *    * *   *  **"); cout<<endl;delay(100);
	cprintf("   *-*-*  *-*-*  *-*-*      *      *    *-*-*     *-*-*    *      *    *   *"); delay(1000);
	clrscr();
	}
}
//**************************************************************************
//			C U R S O R
//**************************************************************************

// SUBJECT MENU
void smcursor()
{  int key;
	while(1)
	{  A:
		while(bioskey(1)==0);
	key=bioskey(0);
	 for(int i=100;i<1000;i+=10)
  {	sound(i);
	delay(1);
	nosound();
	}

	switch(key)
	{
		case 18432 :
						if(smpos==1)
						{ smmenu5();   }
						else if(smpos==2)
						{ smmenu1();  }
						else if(smpos==3)
						{ smmenu2(); }
						else if(smpos==4)
						{ smmenu3();}
						else if(smpos==5)
						{ smmenu4(); }
						else
						goto A;
						break;

		case  20480 :
						if( smpos==1)
						{ smmenu2();  }
						else if(smpos==2)
						{ smmenu3();  }
						else if(smpos==3)
						{ smmenu4();  }
						else if(smpos==4)
						{ smmenu5(); }
						else if(smpos==5)
						{ smmenu1(); }
						else
						goto A;
						break;

		case 7181 : if(smpos==1)
						{  mode='P';
							pcmenu1();
						}
						else if( smpos==2)
						{  mode='C';
							pcmenu1();
						}
						else if( smpos==3)
						matmenu1();
						else if(smpos==4)
						goto A;
						else if(smpos==5)
						{ exit(0); }
						break;
		case 283 : exit(0);                		default :    goto A;

	}
	}
}
//		PHYSICS

void pcursor()
{  int key;
	while(1)
	{  A:
		while(bioskey(1)==0);
	key=bioskey(0);
	 for(int i=100;i<1000;i+=10)
  {	sound(i);
	delay(1);
	nosound();
	}

	switch(key)
	{
		case 18432 :
						if(ppos==1)
						{ pcmenu4();   }
						else if(ppos==2)
						{ pcmenu1();  }
						else if(ppos==3)
						{ pcmenu2(); }
						else if(ppos==4)
						{ pcmenu3();   }
						else
						goto A;
						break;

		case  20480 :
						if( ppos==1)
						{ pcmenu2();  }
						else if(ppos==2)
						{ pcmenu3();  }
						else if(ppos==3)
						{ pcmenu4();  }
						else if(ppos==4)
						{ pcmenu1();  }
						else
						goto A;
						break;

		case 7181 : if(ppos==1)
						 summary();
						else if( ppos==2)
						 auto_mode_phy();
						else if( ppos==3)
						 help();
						else if( ppos==4)
						smmenu1();
						else
						goto A;
						break;
		case 283 : exit(0);
		default :    goto A;

	}
	}
}
// CHEMISTRY
void ccursor()
{  int key;
	while(1)
	{  A:
		while(bioskey(1)==0);
	key=bioskey(0);
	 for(int i=100;i<1000;i+=10)
  {	sound(i);
	delay(1);
	nosound();
	}

	switch(key)
	{
		case 18432 :
						if(cpos==1)
						{ pcmenu4();   }
						else if(cpos==2)
						{ pcmenu1();  }
						else if(cpos==3)
						{ pcmenu2(); }
						else if(cpos==4)
						{ pcmenu3();   }
						else
						goto A;
						break;

		case  20480 :
						if( cpos==1)
						{ pcmenu2();  }
						else if(cpos==2)
						{ pcmenu3();  }
						else if(cpos==3)
						{ pcmenu4();  }
						else if(cpos==4)
						{ pcmenu1();  }
						else
						goto A;
						break;

		case 7181 : if(cpos==1)
						goto A;
						else if( cpos==2)
						 auto_mode_chem();
						else if( cpos==3)
						goto A;
						else if( cpos==4)
						smmenu1();
						else
						goto A;
						break;
		case 283 : exit(0);
		default :    goto A;

	}
	}
}

// MATRIX
void matcursor()
{  int key;
	while(1)
	{  B:
		while(bioskey(1)==0);
	key=bioskey(0);
	 for(int i=100;i<1000;i+=10)
  {	sound(i);
	delay(1);
	nosound();
	}

	switch(key)
	{
		case 18432 :
						if(matpos==1)
						{ matmenu4();   }
						else if(matpos==2)
						{ matmenu1();  }
						else if(matpos==3)
						{ matmenu2(); }
						else if(matpos==4)
						{ matmenu3();   }
						else
						goto B;
						break;

		case  20480 :
						if( matpos==1)
						{ matmenu2();  }
						else if(matpos==2)
						{ matmenu3();  }
						else if(matpos==3)
						{ matmenu4();  }
						else if(matpos==4)
						{ matmenu1();  }
						else
						goto B;
						break;

		case 7181 : if(matpos==1)
						about_matrix();
						else if( matpos==2)
						operation();
						else if( matpos==3)
						inmatmenu1();
						else if( matpos==4)
						smmenu3();
						else
						goto B;
		case 283 : exit(0);
		default :    goto B;

	}
	}
}
//  INMATRIX
void inmatcursor()
{  int key;
	while(1)
	{  A:
		while(bioskey(1)==0);
		key=bioskey(0);
		tone();

	switch(key)
	{
		case 18432 :
						if(inmatpos==1)
						{ inmatmenu6();   }
						else if(inmatpos==2)
						{ inmatmenu1();  }
						else if(inmatpos==3)
						{ inmatmenu2(); }
						else if(inmatpos==4)
						{ inmatmenu3();   }
						else if(inmatpos==5)
						{ inmatmenu4();  }
						else if(inmatpos==6)
						{ inmatmenu5(); }
						else
						goto A;


		case  20480 :
						if( inmatpos==1)
						{ inmatmenu2();  }
						else if(inmatpos==2)
						{ inmatmenu3();  }
						else if(inmatpos==3)
						{ inmatmenu4();  }
						else if(inmatpos==4)
						{ inmatmenu5();  }
						else if(inmatpos==5)
						{ inmatmenu6();  }
						else if(inmatpos==6)
						{ inmatmenu1(); }
						else
						goto A;

		case 7181 : if(inmatpos==1)
						mat_sum();
						else if( inmatpos==2)
						mat_sub();
						else if( inmatpos==3)
						mat_multi();
						else if( inmatpos==4)
						mat_tran();
						else if( inmatpos==5)
						mat_inv();
						else if(inmatpos==6)
						matmenu3();
		default :    goto A;

	}
	}
}


//**************************************************************************
//			  M E N U
//**************************************************************************

//		SUBJECT MENU
void smmenu2()
{  clrscr();
	smpos=2;
	for(int i=0;i<80;i++)
	{
	  textcolor(GREEN);
	  gotoxy(i,0);
	  cprintf("#");
	}        cout<<"\n \n";
	textcolor(CYAN);
	cprintf("*-*-* *-*-* *     *-*-* *-*-* *-*-*   *-*-* *   * *-*-* *-*-* *-*-* *-*-* *-*-*");cout<<endl;
	cprintf("*     *     *     *     *       *     *     *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*   *     *-*   *       *     *-*-* *   * *--*    *   *-*   *       *");cout<<endl;
	cprintf("    * *     *     *     *       *         * *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*-* *-*-* *-*-* *-*-*   *     *-*-* *-*-* *-*-* *-*   *-*-* *-*-*   *");cout<<endl<<"\n \n";
	for(int k=0;k<80;k++)
	{
	  textcolor(GREEN);
	  gotoxy(k,11);
	  cprintf("#");
	}

	cout<<"  \n \n \n \n  \t \t \t";
	textcolor(WHITE);
	textbackground(BLACK);
	cprintf(" 1. P H Y S I C S ");cout<<endl<<"\n  \t \t \t";
	textcolor(BLACK);
	textbackground(GREEN);
	cprintf(" 2. C H E M I S T R Y ");cout<<endl<<"\n  \t \t \t";
	textcolor(WHITE);
	textbackground(BLACK);
	cprintf(" 3. M A T H E M A T I C S ");cout<<endl<<"\n \t \t  \t";
	cprintf(" 4. N E E D  A N Y  H E L P");cout<<endl<<"\n \t \t  \t";
	cprintf(" 5. E X I T");
	smcursor();
}
void smmenu1()
 { clrscr();
	smpos=1;
	for(int i=0;i<80;i++)
	{
	  textcolor(GREEN);
	  gotoxy(i,0);
	  cprintf("#");
	}        cout<<"\n \n";
	textcolor(CYAN);
	cprintf("*-*-* *-*-* *     *-*-* *-*-* *-*-*   *-*-* *   * *-*-* *-*-* *-*-* *-*-* *-*-*");cout<<endl;
	cprintf("*     *     *     *     *       *     *     *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*   *     *-*   *       *     *-*-* *   * *--*    *   *-*   *       *");cout<<endl;
	cprintf("    * *     *     *     *       *         * *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*-* *-*-* *-*-* *-*-*   *     *-*-* *-*-* *-*-* *-*   *-*-* *-*-*   *");cout<<endl<<"\n \n";
	for(int k=0;k<80;k++)
	{
	  textcolor(GREEN);
	  gotoxy(k,11);
	  cprintf("#");
	}

	cout<<"  \n \n \n \n \t \t \t";
	textcolor(BLACK);
	textbackground(GREEN);
	cprintf(" 1. P H Y S I C S ");cout<<endl<<"\n  \t \t \t";
	textcolor(WHITE);
	textbackground(BLACK);
	cprintf(" 2. C H E M I S T R Y ");cout<<endl<<"\n  \t \t \t";
	cprintf(" 3. M A T H E M A T I C S");cout<<endl<<"\n \t  \t \t";
	cprintf(" 4. N E E D  A N Y  H E L P");cout<<endl<<"\n \t \t  \t";
	cprintf(" 5. E X I T");
	smcursor();
}

void smmenu3()
{  clrscr();
	smpos=3;
	for(int i=0;i<80;i++)
	{
	  textcolor(GREEN);
	  gotoxy(i,0);
	  cprintf("#");
	}        cout<<"\n \n";
	textcolor(CYAN);
	cprintf("*-*-* *-*-* *     *-*-* *-*-* *-*-*   *-*-* *   * *-*-* *-*-* *-*-* *-*-* *-*-*");cout<<endl;
	cprintf("*     *     *     *     *       *     *     *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*   *     *-*   *       *     *-*-* *   * *--*    *   *-*   *       *");cout<<endl;
	cprintf("    * *     *     *     *       *         * *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*-* *-*-* *-*-* *-*-*   *     *-*-* *-*-* *-*-* *-*   *-*-* *-*-*   *");cout<<endl<<"\n \n";
	for(int k=0;k<80;k++)
	{
	  textcolor(GREEN);
	  gotoxy(k,11);
	  cprintf("#");
	}

	cout<<"  \n \n \n \n \t \t \t";
	textcolor(WHITE);
	textbackground(BLACK);
	cprintf(" 1. P H Y S I C S ");cout<<endl<<"\n \t  \t \t";
	cprintf(" 2. C H E M I S T R Y ");cout<<endl<<"\n  \t \t \t";
	textcolor(BLACK);
	textbackground(GREEN);
	cprintf(" 3. M A T H E M A T I C S ");cout<<endl<<"\n  \t \t \t";
	textcolor(WHITE);
	textbackground(BLACK);
	cprintf(" 4. N E E D  A N Y  H E L P");cout<<endl<<"\n \t \t  \t";
	cprintf(" 5. E X I T");
	smcursor();
}


void smmenu4()
{  clrscr();
	smpos=4;
	for(int i=0;i<80;i++)
	{
	  textcolor(GREEN);
	  gotoxy(i,0);
	  cprintf("#");
	}        cout<<"\n \n";
	textcolor(CYAN);
	cprintf("*-*-* *-*-* *     *-*-* *-*-* *-*-*   *-*-* *   * *-*-* *-*-* *-*-* *-*-* *-*-*");cout<<endl;
	cprintf("*     *     *     *     *       *     *     *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*   *     *-*   *       *     *-*-* *   * *--*    *   *-*   *       *");cout<<endl;
	cprintf("    * *     *     *     *       *         * *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*-* *-*-* *-*-* *-*-*   *     *-*-* *-*-* *-*-* *-*   *-*-* *-*-*   *");cout<<endl<<"\n \n";
	for(int k=0;k<80;k++)
	{
	  textcolor(GREEN);
	  gotoxy(k,11);
	  cprintf("#");
	}

	cout<<"  \n \n \n \n \t \t  \t";
	textcolor(WHITE);
	textbackground(BLACK);
	cprintf(" 1. P H Y S I C S ");cout<<endl<<"\n \t  \t \t";
	cprintf(" 2. C H E M I S T R Y ");cout<<endl<<"\n  \t \t \t";
	cprintf(" 3. M A T H E M A T I C S ");cout<<endl<<"\n \t \t  \t";
	textcolor(BLACK);
	textbackground(GREEN);
	cprintf(" 4. N E E D  A N Y  H E L P");cout<<endl<<"\n \t \t  \t";
	textcolor(WHITE);
	textbackground(BLACK);
	cprintf(" 5. E X I T");
	smcursor();
}

void smmenu5()
{  clrscr();
	smpos=5;
	for(int i=0;i<80;i++)
	{
	  textcolor(GREEN);
	  gotoxy(i,0);
	  cprintf("#");
	}        cout<<"\n \n";
	textcolor(CYAN);
	cprintf("*-*-* *-*-* *     *-*-* *-*-* *-*-*   *-*-* *   * *-*-* *-*-* *-*-* *-*-* *-*-*");cout<<endl;
	cprintf("*     *     *     *     *       *     *     *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*   *     *-*   *       *     *-*-* *   * *--*    *   *-*   *       *");cout<<endl;
	cprintf("    * *     *     *     *       *         * *   * *   *   *   *     *       *");cout<<endl;
	cprintf("*-*-* *-*-* *-*-* *-*-* *-*-*   *     *-*-* *-*-* *-*-* *-*   *-*-* *-*-*   *");cout<<endl<<"\n \n";
	for(int k=0;k<80;k++)
	{
	  textcolor(GREEN);
	  gotoxy(k,11);
	  cprintf("#");
	}

	cout<<"  \n \n \n \n \t \t \t";
	textcolor(WHITE);
	textbackground(BLACK);
	cprintf(" 1. P H Y S I C S ");cout<<endl<<"\n \t  \t \t";
	cprintf(" 2. C H E M I S T R Y ");cout<<endl<<"\n  \t \t \t";
	cprintf(" 3. M A T H E M A T I C S ");cout<<endl<<"\n \t  \t \t";
	cprintf(" 4. N E E D  A N Y  H E L P");cout<<endl<<"\n \t \t  \t";
	textcolor(BLACK);
	textbackground(GREEN);
	cprintf(" 5. E X I T");
	textcolor(WHITE);
	textbackground(BLACK);
	smcursor();
}
//	PHYSICS/ CHEMISTRY
  void pcmenu1()
{ clrscr(); ppos=1; cpos=1;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n \n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 1. S U M M A R Y  O F  T H E  U N I T");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 2. P R O B L E M   S O L V E R");
  cout<<"\n \n\t\t";
  cprintf(" 3. H E L P  O N  H O W  T U T  C A N  H E L P  Y O U");
  cout<<"\n \n\t\t";
  cprintf(" 4. B A C K  T O  S U B J E C T  L I S T ");
  if(mode=='P')
  pcursor();
  else if(mode=='C')
  ccursor();

}
void pcmenu2()
{ clrscr();    ppos=2;    cpos=2;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. S U M M A R Y  O F  T H E  U N I T");
  cout<<"\n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 2. P R O B L E M   S O L V E R");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 3. H E L P  O N  H O W  T U T  C A N  H E L P  Y O U");
  cout<<"\n \n\t\t";
  cprintf(" 4. B A C K  T O  S U B J E C T  L I S T ");
  if(mode=='P')
  pcursor();
  else if(mode=='C')
  ccursor();

}
 void pcmenu3()
{ clrscr();       ppos=3;      cpos=3;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. S U M M A R Y  O F  T H E  U N I T");
  cout<<"\n \n\t\t";
  cprintf(" 2. P R O B L E M   S O L V E R");
  cout<<"\n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 3. H E L P  O N  H O W  T U T  C A N  H E L P  Y O U");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 4. B A C K  T O  S U B J E C T  L I S T ");
  if(mode=='P')
  pcursor();
  else if(mode=='C')
  ccursor();

}
void pcmenu4()
{ clrscr();         ppos=4;   cpos=4;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. S U M M A R Y  O F  T H E  U N I T");
  cout<<"\n \n\t\t";
  cprintf(" 2. P R O B L E M   S O L V E R");
  cout<<"\n \n\t\t";
  cprintf(" 3. H E L P  O N  H O W  T U T  C A N  H E L P  Y O U");
  textcolor(BLACK);
  textbackground(GREEN);
  cout<<"\n \n \t \t";
  cprintf(" 4. B A C K  T O  S U B J E C T  L I S T ");
  textbackground(BLACK);
  if(mode=='P')
  pcursor();
  else if(mode=='C')
  ccursor();

}

//		MATRIX
 void matmenu1()
{  _setcursortype(_NOCURSOR);
	clrscr(); matpos=1;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n \n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 1. A B O U T  M A T R I C E S");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 2. O P E R A T I O N  O N  M A T R I C E S");
  cout<<"\n \n\t\t";
  cprintf(" 3. P R O B L E M  O N  M A T R I C E S");
  cout<<"\n \n\t\t";
  cprintf(" 4. B A C K  T O  S U B J E C T  L I S T ");
  matcursor();

}
void matmenu2()
{   _setcursortype(_NOCURSOR);
	clrscr();    matpos=2;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. A B O U T  M A T R I C E S");
  cout<<"\n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 2. O P E R A T I O N  O N  M A T R I C E S");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 3. P R O B L E M  O N  M A T R I C E S");
  cout<<"\n \n\t\t";
  cprintf(" 4. B A C K  T O  S U B J E C T  L I S T");
  matcursor();

}
 void matmenu3()
{   _setcursortype(_NOCURSOR);
	clrscr();       matpos=3;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. A B O U T  M A T R I C E S");
  cout<<"\n \n\t\t";
  cprintf(" 2. O P E R A T I O N  O N  M A T R I C E S");
  cout<<"\n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 3. P R O B L E M  O N  M A T R I C E S");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 4. B A C K  T O  S U B J E C T  L I S T");
  matcursor();

}
void matmenu4()
{   _setcursortype(_NOCURSOR);
	clrscr();         matpos=4;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. A B O U T  M A T R I C E S");
  cout<<"\n \n\t\t";
  cprintf(" 2. O P E R A T I O N  O N  M A T R I C E S");
  cout<<"\n \n\t\t";
  cprintf(" 3. P R O B L E M  O N  M A T R I C E S");
  textcolor(BLACK);
  textbackground(GREEN);
  cout<<"\n \n \t \t";
  cprintf(" 4. B A C K  T O  S U B J E C T  L I S T");
  textbackground(BLACK);
  matcursor();

}
//  INMATRIX
void inmatmenu1()
{ clrscr(); inmatpos=1;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n \n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 1. A D D I T I O N  O F  M A T R I X");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 2. S U B T R A C T I O N  O F M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 3. M U L T I P L I C A T I O N  O F  M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 4. T R A N S P O S E  O F  M T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 5. I N V E R S E  O F  M A T R I X");
  cout<<"\n \n \t\t";
  cprintf(" 6. B A C K  T O  M A T R I X  M E N U");
  inmatcursor();

}
void inmatmenu2()
{ clrscr();    inmatpos=2;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. A D D I T I O N  O F  M A T R I X");
  cout<<"\n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 2. S U B T R A C T I O N  O F M A T R I X");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 3. M U L T I P L I C A T I O N  O F  M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 4. T R A N S P O S E  O F  M T R I X");
  cout<<"\n \n \t \t";
  cprintf(" 5. I N V E R S E  O F  M A T R I X");
  cout<<"\n \n \t\t";
  cprintf(" 6. B A C K  T O  M A T R I X  M E N U");
  inmatcursor();

}
 void inmatmenu3()
{ clrscr();       inmatpos=3;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. A D D I T I O N  O F  M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 2. S U B T R A C T I O N  O F M A T R I X");
  cout<<"\n \n \t \t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 3. M U L T I P L I C A T I O N  O F  M A T R I X");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 4. T R A N S P O S E  O F  M T R I X");
  cout<<"\n \n \t \t";
  cprintf(" 5. I N V E R S E  O F  M A T R I X");
  cout<<"\n \n \t\t";
  cprintf(" 6. B A C K  T O  M A T R I X  M E N U");
  inmatcursor();

}
void inmatmenu4()
{ clrscr();         inmatpos=4;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n ";
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n\t\t";
  cprintf(" 1. A D D I T I O N  O F  M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 2. S U B T R A C T I O N  O F M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 3. M U L T I P L I C A T I O N  O F  M A T R I X");
  textcolor(BLACK);
  textbackground(GREEN);
  cout<<"\n \n \t \t";
  cprintf(" 4. T R A N S P O S E  O F  M T R I X");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n \t \t";
  cprintf(" 5. I N V E R S E  O F  M A T R I X");
  cout<<"\n \n \t\t";
  cprintf(" 6. B A C K  T O  M A T R I X  M E N U");
  inmatcursor();

}

 void inmatmenu5()
{ clrscr(); inmatpos=5;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n \n \n \t \t";
  textcolor(WHITE);
  textbackground(BLACK);
  cprintf(" 1. A D D I T I O N  O F  M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 2. S U B T R A C T I O N  O F M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 3. M U L T I P L I C A T I O N  O F  M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 4. T R A N S P O S E  O F  M T R I X");
  cout<<"\n \n\t\t";
  textcolor(BLACK);
  textbackground(GREEN);
  cprintf(" 5. I N V E R S E  O F  M A T R I X");
  textcolor(WHITE);
  textbackground(BLACK);
  cout<<"\n \n \t\t";
  cprintf(" 6. B A C K  T O  M A T R I X  M E N U");
  inmatcursor();

}


 void inmatmenu6()
{ clrscr(); inmatpos=6;
  for(int i=0;i<80;i++)
  {
  gotoxy(i,0);
  textcolor(GREEN+YELLOW+LIGHTRED);
  cprintf("#");
  }
  textcolor(WHITE);
  cout<<"\n";
  cprintf("        #    #   /@\\   * * *  &   &     #    #  * * *  &   &  @   @");
  cout<<endl;
  cprintf("        # \\/ #  @_ _@    *    & \\ &     # \\/ #  * *    & \\ &  @   @");
  cout<<endl;
  cprintf("        #    #  @   @  * * *  &   &     #    #  * * *  &   &  @ @ @");
	cout<<endl;
	for(int j=1;j<=80;j++)
  { gotoxy(j,7);
  textcolor(LIGHTGREEN);
  cprintf("#");
  }
  cout<<"\n \n \n \n \t \t";
  textcolor(WHITE);
  textbackground(BLACK);
  cprintf(" 1. A D D I T I O N  O F  M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 2. S U B T R A C T I O N  O F M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 3. M U L T I P L I C A T I O N  O F  M A T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 4. T R A N S P O S E  O F  M T R I X");
  cout<<"\n \n\t\t";
  cprintf(" 5. I N V E R S E  O F  M A T R I X");
  textcolor(BLACK);
  textbackground(GREEN);
  cout<<"\n \n \t\t";
  cprintf(" 6. B A C K  T O  M A T R I X  M E N U");
  textcolor(WHITE);
  textbackground(BLACK);
  inmatcursor();

}
//***************************************************************************
//	FUNCTION NAME : aboutmatrix
//	DETAILS   	  : IT GIVES DETAILS ABOUT MATRICES
//**************************************************************************
void about_matrix()
{   clrscr();
	char ch;
	ifstream fin("abtmat.txt");
	if(!fin)
	{	cout<<"\n abtmat.txt not found";
		getch();
	}
	while(fin)
	{	fin.get(ch);
		cout<<ch;
	}
	getch();
	matmenu1();
}
//**************************************************************************
//	FUNCTION NAME	: opreation
// DETAILS			: IT
//***************************************************************************
void operation()
{	int n=1,key; char ch;
	ifstream fin("addmat.txt");
	clrscr();
	while(!fin.eof())
	{	fin.get(ch);
		cout<<ch;
	}
	fin.close();
	A:
	while(bioskey(1)==0);
	key=bioskey(0);
	if(key==19200)
	{	if(n!=1) n--;
		else goto A;
	}
	else if(key==19712)
	{	if(n!=4) n++;
		else goto A;
	}
	else if(key==283) matmenu2();
	else goto A;

	clrscr();
	if(n==1){   fin.open("addmat.txt");
					fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==2){ fin.open("submat.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	 }
	else if(n==3){ fin.open("multimat.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==4){ fin.open("transmat.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	goto A;
}
//**************************************************************************
//	FUNCTION NAME : mat_sum
// DETAILS		  : IT
//**************************************************************************
void mat_sum()
{  _setcursortype(_NORMALCURSOR);
	clrscr();
	int row,col,no,ar[100][10][10],sum[10][10];   int i,j,k;
	cout<<"\n enter number of matix :- ";
	cin>>no;
	cout<<"\n enter row and column of "<<no<<" matrix :- ";
	cin>>row>>col;

	for(k=0;k<no;k++)
	{	cout<<"\n enter value of matrix-"<<k+1<<" :- ";
		for(i=0;i<row;i++)
		{	for(j=0;j<col;j++)
				cin>>ar[k][i][j];
		}
	}
	for(i=0;i<row;i++)
	{	for(j=0;j<col;j++)
			sum[i][j]=0;
	}

	for(k=0;k<no;k++)
	{	for(i=0;i<row;i++)
		{  for(j=0;j<col;j++)
				sum[i][j]+=ar[k][i][j];
		}
	}
	cout<<" sum of two matrix is...."<<" ";

	for(i=0;i<row;i++)
	{  cout<<"\n";
		for(j=0;j<col;j++)
		{
			cout<<"  "<<sum[i][j];

		}

	}
	getch();
	_setcursortype(_NOCURSOR);
	inmatmenu1();
}
//**************************************************************************
//	FUNCTION NAME : mat_sub
//	DETAILS		  : IT
//***************************************************************************
void mat_sub()
{  _setcursortype(_NORMALCURSOR);
 clrscr();
 int row,col,arr1[10][10],arr2[10][10],sub[10][10];
 int i,j,k;


 cout<<"\n enter row and column of matrix  :- ";
 cin>>row>>col;

 cout<<"\n enter value of matrix-1 :- ";
 for(i=0;i<row;i++)
 {  for(j=0;j<col;j++)
	 {
		cin>>arr1[i][j];
	 }
 }

 cout<<"\n enter value of matrix-2 :- ";

 for(i=0;i<row;i++)
 {  for(j=0;j<col;j++)
	 {
		cin>>arr2[i][j];
	 }
 }
 cout<<"\n value of matrix-1 :- ";

 for(i=0;i<row;i++)
 {  cout<<"\n";
	for(j=0;j<col;j++)
	{
	  cout<<"  "<<arr1[i][j];

	}
 }

 cout<<"\n value of matrix-2 :- ";

 for(i=0;i<row;i++)
 {  cout<<"\n";
	for(j=0;j<col;j++)
	{
	  cout<<"  "<<arr2[i][j];

	}
 }

 for(i=0;i<row;i++)
 {
	for(j=0;j<col;j++)
	{
		sub[i][j]=0;
	}
 }
 cout<<"\n subtraction value is :-  ";
 for(i=0;i<row;i++)
 {
	for(j=0;j<col;j++)
	{
	  sub[i][j]=arr1[i][j]-arr2[i][j];
	}
 }

 cout<<"\n";
 for(i=0;i<row;i++)
 { cout<<"\n";
	for(j=0;j<col;j++)
	{
	  cout<<"  "<<sub[i][j];

	}
 }
 getch();
 _setcursortype(_NOCURSOR);
 inmatmenu2();

}
//***************************************************************************
//	FUNCTION NAME : mat_multi
//	DETAILS    	  : IT
//***************************************************************************
void mat_multi()
{  _setcursortype(_NORMALCURSOR);
	int i,j,k;
  clrscr();
 int row1,col1, row2,col2,arr1[10][10],arr2[10][10],multi[10][10];

 cout<<"\n enter row and column of matrix 1st :- ";
 cin>>row1;cout<<"\t\t\t\t\t\b";cin>>col1;

 cout<<"\n enter row and column of matrix 2nd :- ";
 cin>>row2;  cout<<"\t\t\t\t\t\b";cin>>col2;
 if(col1==row2)
 {
	cout<<"\n enter value of matrix 1st ("<<row1<<"x"<<col1<<") :-  ";

	for(i=0;i<row1;i++)
  {
	 for(j=0;j<col1;j++)
	 {
		cin>>arr1[i][j];
	 }

  } 	cout<<"\n enter value of matrix 2nd ("<<row2<<"x"<<col2<<") :-  ";

	for(i=0;i<row2;i++)
  {
	 for(j=0;j<col2;j++)
	 {
		cin>>arr2[i][j];
	 }

  }
	clrscr();
	cout<<"\n value of matrix 1st :-  ";

	for(i=0;i<row1;i++)
  {
		cout<<"\n";
		for(j=0;j<col1;j++)
		{
			cout<<"  "<<arr1[i][j];
		}
  }
	cout<<"\n value of matrix 2nd :-  ";

	for(i=0;i<row2;i++)
  {
	  cout<<"\n";

	 for(j=0;j<col2;j++)
	 {
		cout<<"  "<<arr2[i][j];

	 }

  }

  cout<<"\n product of two matrix is :-  ";

	for(i=0;i<row1;i++)
	{
		cout<<"\n";
		for(j=0;j<col2;j++)
		{
			multi[i][j]=0;
			for(k=0;k<col1;k++)
				multi[i][j]=multi[i][j]+arr1[i][k]*arr2[k][j];
				cout<<"  "<<multi[i][j];

		}
	}
  }
	else
	{
		cout<<"\n matrix are not compatile for multiplication ";
	}
	getch();
	_setcursortype(_NOCURSOR);
	inmatmenu3();
 }
//***************************************************************************
//	FUNCTION NAME : mat_tran
//	DETAILS		  : IT
//***************************************************************************
void mat_tran()
{	_setcursortype(_NORMALCURSOR);
 clrscr();char inp[3];
 int row,col,ar[10][10],tra[10][10];   int i,j,k;
 cout<<"\n enter row and column of matrix :- ";
 cin>>row;cout<<"\t\t\t\t\t\b\b\b\b";cin>>col;

 cout<<"\n enter value of matrix ("<<row<<"x"<<col<<") :-  ";

  for(i=0;i<row;i++)
  {
	 for(j=0;j<col;j++)
	 { A:
		cin>>inp;
		if( isdigit(inp[0]) )
		ar[i][j]=atoi(inp);
		else
		{ cout<<"\n Enter a number:";
			goto A;
		}

		tra[i][j]=0;
		tra[i][j]=ar[i][j];
	 }

  }
	clrscr();
	cout<<"\n value of matrix :-  ";

	for(i=0;i<row;i++)
  {
	  cout<<"\n";

	 for(j=0;j<col;j++)
	 cout<<"  "<<ar[i][j];
  }


  cout<<"\n\n transpose of matrix is :-  ";

	for(i=0;i<col;i++)
		{
			cout<<"\n";
			for(j=0;j<row;j++)
			 cout<<"  "<<tra[j][i];
		}
		getch();
		_setcursortype(_NOCURSOR);
		inmatmenu4();

}
//***************************************************************************
//	FUNCTION NAME : mat_inv
//	DETAILS		  : IT
//***************************************************************************
void mat_inv()
{  _setcursortype(_NORMALCURSOR);
	clrscr();
	float a[3][3],b[3][3],c[3][3],det;
	cout<<"\n enter matrix element(3x3)"<<"\n";
	for(int i=0;i<3;i++)
	 for(int j=0;j<3;j++)
	  cin>>a[i][j];

	cout<<"\n value of matrix is :- ";
	 for(i=0;i<3;i++)
	 {  cout<<"\n";
		 for(j=0;j<3;j++)
		 { cout<<"  "<<a[i][j];
		 }
	 }
	b[0][0]=+(a[1][1]*a[2][2]-a[1][2]*a[2][1]);
	b[0][1]=-(((a[2][2]*a[1][0])-(a[2][0]*a[1][2])));
	b[0][2]=a[1][0]*a[2][1]-a[1][1]*a[2][0];
	b[1][0]=-(a[0][1]*a[2][2]-a[0][2]*a[2][1]);
	b[1][1]=a[2][2]*a[0][0]-a[0][2]*a[2][0];
	b[1][2]=-((a[0][0]*a[2][1])-(a[2][0]*a[0][1]));
	b[2][0]=a[0][1]*a[1][2]-a[0][2]*a[1][1];
	b[2][1]=-(a[0][0]*a[1][2]-a[0][2]*a[1][0]);
	b[2][2]=a[0][0]*a[1][1]-a[0][1]*a[1][0];

	for(int q=0;q<3;q++)
	{
		for(int w=0;w<3;w++)
		{
			c[q][w]=b[w][q];
		}
	}
	det=b[0][0]*a[0][0]+b[0][1]*a[0][1]+b[0][2]*a[0][2];
	if(det==0)
	{
		cout<<"\n inverse not possible because determinant of matrix is 0";
		getch();
	  //	return;
	  _setcursortype(_NOCURSOR);
	  inmatmenu5();
	}
	else
	{  cout<<"\n\n the inverse of matrix is:\n";
		for(int s=0;s<3;s++)
		{ cout<<endl;
		 for(int d=0;d<3;d++)
			 {cout.setf(ios::fixed);
			  cout<<setw(8)<<setprecision(2)<<(c[s][d])/det<<" ";
			 }
		}
	}
	getch();
  //	return;
  _setcursortype(_NOCURSOR);
  inmatmenu5();
}

//**************************************************************************
//	FUNCTION NAME : tone
//	DETAILS	      : IT PRODUCES A SOUND ON KEY PRESS OR SCREEN CHANGE
//***************************************************************************
void tone()
{     for(int k=100;k<1000;k+=10)
		{	sound(k);
	delay(1);
	nosound();
		}
}
//***************************************************************************
//	FUNCTION NAME : auto_mode
//	DETAILS		  : IT TAKES THE PROBLEM FROM THE USER AND USES THE
//					    THE OTHER FUNCTIONS TO SOLVE THE PROBLEM
//***************************************************************************
void auto_mode_phy()
{  clrscr();
	_setcursortype(_NORMALCURSOR);
	cout<<"\n ENTER YOUR QUESTION:\n\n\n";
	gets(tmp);
	getwords();
	shape_identifier();
	data_extracter_phy();
	strlwr(tmp);
	getwords();
	aim_phy();
	function_selecter_phy();
	_setcursortype(_NOCURSOR);
	pcmenu2();
}
void auto_mode_chem()
{	clrscr();
	_setcursortype(_NORMALCURSOR);
	cout<<"\n ENTER YOUR QUESTION:\n\n\n";
	gets(tmp);
	getwords();
	pack_identifier();
	data_extracter_chem();
	strlwr(tmp);
	getwords();
	aim_chem();
	_setcursortype(_NOCURSOR);
	pcmenu2();
}
//***************************************************************************
//	MAIN() FUNCTION
//***************************************************************************
void main()
{  clrscr();
   cout.setf(ios::scientific);
  _setcursortype(_NOCURSOR);
  clrscr();
  scr();
  jnv();
  tone();
  delay(2200);
  name();
  tone();
  delay(2200);
  tone();
 // help();
 smmenu1();
 tone();
}