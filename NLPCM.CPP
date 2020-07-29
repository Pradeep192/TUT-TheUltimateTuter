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
double z,M,a,density,re;
//**************************************************************************
//	GLOBAL VARIABLES
//**************************************************************************
char tmp[1001],input[100][31],unit[10],cstring[30];
int itf=0,tw,to_find[5]={-2,-2,-2,-2,-2};
char d=196,sq=253,cpi=227,ceps=238,cs=229,ctau=226,ctheta=233,cphi=237;
char cdell=127,carr=26,sqr=251,line=179,pack[30]="null",mode;
//**************************************************************************
//	FUNCTION PROTOTYPES
//**************************************************************************
void getwords();
void convert_to_SI(quantity&);
void data_initializer_chem(char[]);
void data_extracter_chem();
void case_selecter_chem(int);
int case_selecter_chem(char[]);
void cell_density(char[]);
void pack_identifier();
void radius_edge(char[]);
int spell_check(char[],char[]);
void aim_chem();
int aim_chem(char[]);
void aim_chem(char[],int&,int);
void remove_punctuation(char[]);
void function_selecter_chem();
int sign_matcher_full_chem(signature);
int sign_matcher_inc_chem(signature);
void function_caller_chem(int);
void auto_mode_chem();
void initializer();
void variable_printer();
void tone();

void csummary();
void summary();
void help(int);
void operation();
void aboutmatrix();

void smcursor();
void pcursor();
void ccursor();
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
	re=M=density=a=z=12321;
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

	if(M!=12321) cout<<"M="<<M<<" kg_mol(-1)  ";
	if(a!=12321) cout<<"a="<<a<<" m  ";
	if(density!=12321) cout<<"d="<<density<<" kg_m(-3)  ";
	if(re!=12321) cout<<"r="<<re<<" m ";
	if(z!=12321) cout<<"z="<<z<<" ";
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
	{
		if( cmp[0]=='d' && cmp[1]=='a')
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
		else if( strcmp(cmp,"u")==0 || strcmp(cmp,"g_mol(-1)")==0)
		{	A.getdata(A.retmag()*1e-3,"kg_mol(-1)");
			break;
		}
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
		cout<<"\n\n The density is:"<<density<<" kg_m(-3)"; getch();
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
		cout<<"\n\n The molar mass is:"<<M<<" kg_mol(-1)"; getch();
	}
	else if( strcmp(find,"nature")==0 )
	{  cout<<"     d x a^3 x N_a";
		cout<<"\n\t\t\t  z= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t  M";
		getch();
		cout<<"\n\n\t\t\t     ("<<density<<") x ("<<a<<")^3 x ("<<N_a<<')';
		cout<<"\n\t\t\t  z= "<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d<<d;
		cout<<"\n\t\t\t\t\t  ("<<M<<')';  getch();
		z=density*a*a*a*N_a/M;
		cout<<"\n z="<<z;
		z=floor(z);
		if(z==1)
		cout<<"\n\n The unit cell is in simple packing.";
		else if(z==2)
		cout<<"\n\n The unit cell is in body centred cubic packing.";
		else if(z==4)
		cout<<"\n\n The unit cell is in face centred cubic packing (/ccp/hcp).";
		else {cout<<"\n no match found  "<<z; }
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
//###########################################################################
//###########################################################################
//***************************************************************************
//***************************************************************************
//	FUNCTION NAME : data_extracter_chem
//	DETAILS		  :
//***************************************************************************
void data_extracter_chem()
{	for(int j,i=0;i<tw;++i)
	{	if( isdigit(input[i][0]) || isdigit(input[i][1]) )
		{
			data.getdata( atof(input[i]), input[i+1] );
			cout<<'\n';data.print_quantity();
			convert_to_SI(data);
			data.print_quantity();
			if( strcmp(data.retunit(),"kg_mol(-1)")==0 )
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
}
//***************************************************************************
//	FUNCTION NAME : pack_identifier
//	DETAILS		  :
//***************************************************************************
void pack_identifier()
{	for(int i=0;i<tw;++i)
	{	if( !strcmp(input[i],"fcc") || !strcmp(input[i],"hcp") || !strcmp(input[i],"ccp") )
		{	strcpy(pack,"fcc"); z=4; qsign.sign[3]++; }
		else if( strcmp(input[i],"hexagonal")==0 && strcmp(input[i+1],"close")==0 )
		{	strcpy(pack,"fcc"); z=4; qsign.sign[3]++; }
		else if( !strcmp(input[i],"bcc") )
		{	strcpy(pack,"bcc"); z=2; qsign.sign[3]++;  }
		else if( strcmp(input[i],"cubic")==0 )
		{	if(!strcmp(input[i-1],"face-centred") )
			{	strcpy(pack,"fcc"); z=4; qsign.sign[3]++;  }
			else if( strcmp(input[i+1],"close")==0 )
			{	strcpy(pack,"fcc"); z=4; qsign.sign[3]++;  }
			else if(strcmp(input[i-1],"body-centred")==0 )
			{	strcpy(pack,"bcc"); z=2; qsign.sign[3]++;  }
			else if( !strcmp(input[i-1],"simple") )
			{	strcpy(pack,"simple"); z=1; qsign.sign[3]++;  }
		}
	}
}
//***************************************************************************
//	FUNCTION NAME : aim
//	DETAILS 		  : IT SEARCHES FOR WHAT IS TO BE CALCULATED
//***************************************************************************
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
	words in; f=0;
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
			break;
		}
		else if( strcmp(str,"mass")==0 )
		{	if( !strcmp(input[i-1],"molar") || !strcmp(input[i-1],"molecular") || !strcmp(input[i-1],"atomic") )
			{	to_find[itf++]=0;
				qsign.sign[0]++; f=1; break;
			}
		}
	}
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

void csummary()
{	int n=1,key; char ch;
	ifstream fin("CSUM1.txt");
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
	else if(key==283) pcmenu1();
	else goto A;

	clrscr();
	if(n==1){   fin.open("CSUM1.txt");
					fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==2){ fin.open("CSUM2.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	 }
	else if(n==3){ fin.open("CSUM3.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==4){ fin.open("CSUM4.txt");
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
void help(int num)
{	int n=1,key; char ch;
	ifstream fin("UMAN1.txt");
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
	{	if(n!=3) n++;
		else goto A;
	}
	else if(key==283 || key==7181)
	{  if(num==0) return;
		if(num==1) smmenu3();
	}
	else goto A;

	clrscr();
	if(n==1){   fin.open("UMAN1.txt");
					fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	}
	else if(n==2){ fin.open("UMAN2.txt");
						fin.seekg(0);
					while(!fin.eof())
					{	fin.get(ch);
						cout<<ch;
					}
					fin.close();
	 }
	 else if(n==3){ fin.open("UMAN3.txt");
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
{  randomize();
	clrscr();
	int i;
	while(bioskey(1)==0)
	{  i=random(1000);
		delay(200);
	textcolor(i);
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
	exit(0);
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
						help(1);
						else if(smpos==5)
						{ exitdesign(); }
						break;
		case 283 : exitdesign();                		default :    goto A;

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
						{ pcmenu3();   }
						else if(ppos==2)
						{ pcmenu1();  }
						else if(ppos==3)
						{ pcmenu2(); }
						else
						goto A;
						break;

		case  20480 :
						if( ppos==1)
						{ pcmenu2();  }
						else if(ppos==2)
						{ pcmenu3();  }
						else if(ppos==3)
						{ pcmenu1();  }
						else
						goto A;
						break;

		case 7181 : if(ppos==1)
						 summary();
						else if( ppos==2)
						 goto A;
						else if( ppos==3)
						smmenu1();
						else
						goto A;
						break;
		case 283 : exitdesign();
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
						{ pcmenu3();   }
						else if(cpos==2)
						{ pcmenu1();  }
						else if(cpos==3)
						{ pcmenu2(); }
						else
						goto A;
						break;

		case  20480 :
						if( cpos==1)
						{ pcmenu2();  }
						else if(cpos==2)
						{ pcmenu3();  }
						else if(cpos==3)
						{ pcmenu1();  }
						else
						goto A;
						break;

		case 7181 : if(cpos==1)
						csummary();
						else if( cpos==2)
						 auto_mode_chem();
						else if( cpos==3)
						smmenu2();
						else
						goto A;
						break;
		case 283 : exitdesign();
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
		case 283 : exitdesign();
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
						goto A;
						else if( inmatpos==2)
						goto A;
						else if( inmatpos==3)
						goto A;
						else if( inmatpos==4)
						goto A;
						else if( inmatpos==5)
						goto A;
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
  cprintf(" 3. B A C K  T O  S U B J E C T  L I S T ");
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
	cprintf(" 3. B A C K  T O  S U B J E C T  L I S T ");
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
  cprintf(" 3. B A C K  T O  S U B J E C T  L I S T ");
  textcolor(WHITE);
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
  cprintf(" 4. T R A N S P O S E  O F  M A T R I X");
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
  cprintf(" 4. T R A N S P O S E  O F  M A T R I X");
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
  cprintf(" 4. T R A N S P O S E  O F  M A T R I X");
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
  cprintf(" 4. T R A N S P O S E  O F  M A T R I X");
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
  cprintf(" 4. T R A N S P O S E  O F  M A T R I X");
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
  cprintf(" 4. T R A N S P O S E  O F  M A T R I X");
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
void auto_mode_chem()
{	clrscr();
	_setcursortype(_NORMALCURSOR);
	initializer();
	cout<<"\n ENTER YOUR QUESTION:\n\n\n";
	gets(tmp);
	getwords();
	data_extracter_chem();
	pack_identifier();
	strlwr(tmp);
	getwords();
	aim_chem();
	variable_printer();
	function_selecter_chem();
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
  help(0);
 smmenu1();
 tone();
 getch();
}