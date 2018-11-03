
import java.io.*;
import java.util.Arrays;

public class ExtractCDSs{

public static void main(String[] args){

try
{
	ExtractCDSs FAobg=new ExtractCDSs();
	
	String WorkDir="/N/dc2/scratch/xw63/Maps/";
	
	String MapFile = args[0];
	String OutputDir = args[1];

	BufferedReader brMap=new BufferedReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(new File(MapFile))),"utf-8"),5*1024*1024); // reading file using 5M buffer;
		
	int Gn=200000;
	String[] b = new String[100];	
	String [] Reference=new String[Gn];
	String [] Major=new String[Gn];
	String [] Minor=new String[Gn];
	String [] ps = new String[Gn]; //Minor Allele freq 
	String [] qs = new String[Gn]; //Major Allele freq 
	String [] Strand = new String[Gn]; //Gene strand 
	int [] Start = new int[Gn]; //Gene start 
	int [] End = new int[Gn]; //Gene end 

	for (int g=1;g<Gn;g++) //Is a coding sequence
	{
			Reference[g]="";
			Major[g]="";
			Minor[g]="";
			ps[g]="";
			qs[g]="";
	}
	System.out.println("Reading annotated map files:\n");
	int i=0;
	String lineMap="";
	
	while((lineMap=brMap.readLine()) != null)
	{		
		b = lineMap.replace(" ", "\t").split("\t");
		String GeneType=b[1]; //Gene type
		i++;
		if (GeneType.equals("CDS")) //skip the first/last/empty line
		{
			int g=Integer.parseInt(b[0].replace("Gene","")); //Gene No.
			String Type=b[1];//Type of gene
			Strand[g]=b[2];//Strand of DNA
			String Ref=b[3];
			String RefChk=b[4];
			int S=Integer.parseInt(b[5]); //Scaffold
			int L= Integer.parseInt(b[6]); //Position			
			String Ref1=b[7];
			String Maj=b[8];
			String Min=b[9];
			
	//Gene    Type    Strand  Ref0    RefChk  Sca     Pos     Ref     Maj     Min     CovStr  MJ_FREQ VR_FREQ ERROR   NULL_ER NULL_ER2        F_STAT  MM_FREQ Mm_FREQ mm_FREQ HETERO  PLR     HWE_LR  GOF     EF_CHRM IND_INC IND_CUT REF_BIAS   P_REF_BS        BEST_LL group   P       HLR     pi      H       Fis
	
			//Minor allele frequencies (MAF)
			if(b[31].equals(".")){b[31]="0.0";	}	
			double MAF = Double.parseDouble(b[31]); //MAF:	
			String p=String.valueOf(Math.round((1-MAF)*10000)/10000.0);
			String q=String.valueOf(Math.round(MAF*10000)/10000.0);
					
			if (Ref.equals("N")) {Ref=Maj;}
			if (Ref.equals("N")) {Ref=Min;}
			if (Maj.equals("N")) {Maj=Ref;}
			if (Min.equals("N")) {Min=Ref;}
				
			if(Strand[g].equals("+"))
			{
				Reference[g]=Reference[g]+" "+Ref;
				Major[g]=Major[g]+" "+Maj;
				Minor[g]=Minor[g]+" "+Min;
				ps[g]=ps[g]+" "+p;
				qs[g]=qs[g]+" "+q;
			}
			else
			{
				int N1="TAGCN".indexOf(Ref);
				int N2="TAGCN".indexOf(Maj);
				int N3="TAGCN".indexOf(Min);
				Ref="ATGCN".substring(N1, N1+1); //compliment
				Maj="ATGCN".substring(N2, N2+1); //compliment
				Min="ATGCN".substring(N3, N3+1); //compliment
				
				Reference[g]=Ref+" "+Reference[g]; //Reverse
				Major[g]=Maj+" "+Major[g];//Reverse
				Minor[g]=Min+" "+Minor[g];//Reverse
				ps[g]=p+" "+ps[g];//Reverse
				qs[g]=q+" "+qs[g];//Reverse
			}
			if(Reference[g].length()==1){Start[g]=L;}
		}
 	}

	System.out.println("Writing CDS files in:"+OutputDir+"\n");

	int Gn1=0;
	for (int g=1;g<Gn;g++) //Is a coding sequence
	{
		if (!Reference[g].isEmpty())
		{
			Gn1++;
			String CDSName="Gene"+g;
			String CDSFileName=OutputDir+CDSName+".txt5";
			FileWriter fw3=new FileWriter(CDSFileName);
			BufferedWriter bw3=new BufferedWriter(fw3);			

			System.out.println(CDSName+"\n");
			
			End[g]=Start[g]+Reference[g].length()-1;
			
			bw3.write(">"+CDSName+"-"+Start[g]+"-"+End[g]+"-("+Strand[g]+")\n"); 
			bw3.write(Reference[g]+"\n"); 	
			bw3.write(Major[g]+"\n"); 	
			bw3.write(Minor[g]+"\n"); 	
			bw3.write(ps[g]+"\n"); 	
			bw3.write(qs[g]+"\n\n");
			bw3.close();
		}
	}	
	System.out.println(Gn1+" CDS files were saved in:"+OutputDir+"\n");
}
catch(Exception e)
{
	System.out.println("\nThis program is to ExtractCDSs\n"); 
	System.out.println("\nUsage: Java -cp ./ ExtractCDSs\n"); 
	System.out.println(e);
	e.printStackTrace();
}
}
public static int countOf (String s, String c) {
    return s.length() - s.replace(c, "").length();
}
}
/*
======================================================================
ExtractCDSs

Usage: ExtractCDSs using Gene Frequency Data and gene annotation file:
 
Java -cp ./ ExtractCDSs <Map-GeneFreqData> <GFF> > output
=====================================================================
Written by:                   
Xiaolong Wang @ Ocean University of China 
email: xiaolong@ouc.edu.cn
website: http://www.DNAplusPro.com
=====================================================================
In hope useful in genomics and bioinformatics studies.
This software is released under GNU/GPL license
Copyright (c) 2018
1. Lynch Lab, CME, Biodesign, Arizona State University
2. Lab of Molecular and Computational Biology, Ocean University of China

You may use it or change it as you wish. 
We will be glad to know if you find it is useful. 
=====================================================================
*/
