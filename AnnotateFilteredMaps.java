
import java.io.*;
import java.util.Arrays;

public class AnnotateFilteredMaps{

public static void main(String[] args)
{
	int Tn=5; //Number of years
	int i=0;
	int j=0;
	int n=0; 
	int S=0; //Scaffold Number
	int L=0; //Locus Position	
	int TotalSites=0;
	String OStr="";
	String [] a = new String [50];
	String lineMap ="";
	
	int Sn=498; //Maximum Number of scaffolds
	int [] Ln = new int [Sn];//Maximum Number of nucleotides in a scaffold
	int Gn=2000000; //Maximum Number of genes
	char [][] RefGenome=new char [Sn][];
	int [][] AnoGene=new int [Sn][];
	int [][] AnoType=new int [Sn][];
	String [][] AnoStrand = new String [Sn][] ; //Gene strand: +: sense; -: antisense

try
{	
	AnnotateFilteredMaps FAobj=new AnnotateFilteredMaps();
	
	String WD= args[0];
	int L0=WD.length();
	if (WD.substring(L0).equals("/")){WD=WD.substring(0,L0-1);}
	
	
	String AnnotatedMapFile = WD+"/"+args[1];	
	BufferedReader brRef = new BufferedReader(new InputStreamReader( new BufferedInputStream(new FileInputStream(new File(AnnotatedMapFile))),"utf-8"),5*1024*1024);// reading file using 5M buffer 
	
	String FilteredMapFile = WD+"/"+args[2];	
	BufferedReader brFilMap = new BufferedReader(new InputStreamReader( new BufferedInputStream(new FileInputStream(new File(FilteredMapFile))),"utf-8"),5*1024*1024);// reading file using 5M buffer 
	
	String HeadFile = WD+"/"+args[3];	
	BufferedReader brHead = new BufferedReader(new InputStreamReader( new BufferedInputStream(new FileInputStream(new File(HeadFile))),"utf-8"),16*1024);// reading file using 16k buffer 
	
	System.out.print("\nReading the head file:\n"+HeadFile+"\n");
	System.out.print("Sca\t");
	System.out.print("Length\n");
	String lineHead="";
	while((lineHead=brHead.readLine()) != null)
	{
		a = lineHead.replace(" ", "\t").split("\t");
		
		//skip the first/last/empty line
		if (lineHead.indexOf("scaffold_")>=0) 
		{					
			S=Integer.parseInt(a[0].substring(9)); //Scaffold
			Ln[S]= Integer.parseInt(a[1])+1; //Length of the scaffold
			RefGenome[S]=new char [Ln[S]];
			AnoGene[S]=new int [Ln[S]];
			AnoType[S]=new int [Ln[S]];
			AnoStrand[S] = new String [Ln[S]] ; 
			System.out.print(S+"\t");
			System.out.print(Ln[S]-1+"\n");
		}
	}
	brHead.close();
 		 			 	
 	System.out.print("Reading the Annotated Map:\n"+AnnotatedMapFile+"\n");
	
	System.out.print("Sca\tPos\tRef\tGene\tType\tStrand\n");
	String GeneStr="";
	String TypeStr="";
	String StrandStr="";
	int GeneNo=0;
	int TypeNo=0;
	while((lineMap=brRef.readLine()) != null)
	{
		i++;
		a=lineMap.toUpperCase().trim().replace(" ", "\t").split("\t");
	
		if (i>1&&a.length>2&&(!(a[0]==null||a[0].isEmpty()))) 
		{//skip the first line
			S=Integer.parseInt(a[0]); //Scaffold
			L=Integer.parseInt(a[1]); //Position			
			RefGenome[S][L]=a[2].charAt(0);
			GeneStr="";
			TypeStr="";
			StrandStr="";
			GeneNo=0;
			TypeNo=0;
			if(a.length==6)
			{
				if (!(a[3]==null||a[3].isEmpty())) 
				{ 
					GeneNo=Integer.parseInt(a[3]);
				}
				if (!(a[4]==null||a[4].isEmpty())) 
				{
					if (a[4].equals("CDS")) {TypeNo=1;} 
					if (a[4].equals("exon")) {TypeNo=2;}
					if (a[4].equals("intron")) {TypeNo=3;}
					if (a[4].equals("UTR3")) {TypeNo=4;}
					if (a[4].equals("UTR5")) {TypeNo=5;}
				}
				if (!(a[5]==null||a[5].isEmpty()))
				{ 
					StrandStr=a[5];
				}
			}
			AnoGene[S][L]=GeneNo;
			AnoType[S][L]=TypeNo;
			AnoStrand[S][L]=StrandStr==""?".":StrandStr;
			GeneStr=GeneNo==0?".":"Gene"+String.valueOf(GeneNo);
			TypeStr=TypeNo==0?".":String.valueOf(TypeNo)+"-"+a[4];
			if (i<5||(i%1000000)==0) //print the first 5 lines
			{
				System.out.print(i+" lines: Scaffold-");
				System.out.print(S+"\t");
				System.out.print(L+"\t");
				System.out.print(a[2]+"\t");
				System.out.print(GeneStr+"\t");
				System.out.print(TypeStr+"\t");
				System.out.print(StrandStr+"\n");
			}
		}
	}

//Show the last line of the Annotated Map
	System.out.print("......\n");
	System.out.print(S+"\t");
	System.out.print(L+"\t");
	System.out.print(RefGenome[S][L]+"\t");
	System.out.print(GeneStr+"\t");
	System.out.print(TypeStr+"\t");
	System.out.print(StrandStr+"\n");
	                       			
	System.out.print("Reading the Filtered Map:\n"+FilteredMapFile+"\n");
		
	String O0 = FilteredMapFile+".Annotated.txt";
	
	System.out.print("Write to file:\n"+O0+"\n");	
		
	BufferedWriter bw0 = new BufferedWriter(new OutputStreamWriter( new BufferedOutputStream(new FileOutputStream(new File(O0))),"utf-8"),1024*1024); // Writing file using 1M buffer 
	
	i=0;
	String OutStr="";
	while((lineMap=brFilMap.readLine()) != null)
	{
		if(lineMap==null){break;}
		if(lineMap.trim().isEmpty()){continue;}
		
		i++;
		a=lineMap.toUpperCase().replace(" ", "\t").split("\t");
			
		if (i==1) //the headline
		{
			OutStr="Gene\tType\tStrand\tRef0\tRefChk\t"+lineMap.trim()+"\n";
		}
		else
		{	
			S=Integer.parseInt(a[0]); //Scaffold
			L=Integer.parseInt(a[1]); //Position
			char R0=RefGenome[S][L];
			char R1=a[3].charAt(0);
			String R01=(R1==R0)?"":"*";
			TotalSites++;
			GeneNo=AnoGene[S][L];
			TypeNo=AnoType[S][L];
			StrandStr=AnoStrand[S][L];
			GeneStr=GeneNo==0?".":"Gene"+String.valueOf(GeneNo);
			TypeStr=".";
			if (TypeNo==1) {TypeStr="CDS";} 
			if (TypeNo==2) {TypeStr="exon";}
			if (TypeNo==3) {TypeStr="intron";}
			if (TypeNo==4) {TypeStr="UTR3";}
			if (TypeNo==5) {TypeStr="UTR5";}
			
			OutStr=GeneStr+"\t"+TypeStr+"\t"+StrandStr+"\t"+R0+"\t"+R1+R01+"\t"+lineMap.trim()+"\n";
			
		}
		bw0.write(OutStr);
		if (i<5||(i%1000000)==0) //print the first 3 lines
		{
			System.out.print(i+" lines: ");
			System.out.print(OutStr);
		}
	 }
	System.out.print("......\n");
	System.out.print(OutStr);
	bw0.close();
}
catch(Exception e)
{
	System.out.println("\nThis program is to AnnotateFilteredMaps\n"); 
	System.out.println("\nUsage: Java -cp ./ AnnotateFilteredMaps <workdir> <FilteredMapFile> <AnnotatedMapFile> <HeadFile>\n"); 
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
AnnotateFilteredMaps

Usage: AnnotateFilteredMaps using Gene Frequency Data and gene annotation file:
 
Java -cp ./ AnnotateFilteredMaps <Map-GeneFreqData> <GFF> > output
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
