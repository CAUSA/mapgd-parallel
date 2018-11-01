
import java.io.*;
import java.util.Arrays;

public class CombineRefGFF{

public static void main(String[] args){

try
{
	CombineRefGFF FAobj=new CombineRefGFF();
	
	String WD="/N/dc2/scratch/xw63/Maps-0.4.34/";
	String HeadFile = WD+"/PA2013.map.idx";	
	String RefGenomeFile = WD+"PA42.4.1.fasta";	
	String GFF = WD+"PA42.4.1.gff";
	String AnoMap = WD+"PA42.4.1-Annotated.map";
	
	BufferedReader brRefGenome = new BufferedReader(new InputStreamReader( new BufferedInputStream(new FileInputStream(new File(RefGenomeFile))),"utf-8"),5*1024*1024);// reading file using 5M buffer 
 		
	BufferedReader brGFF = new BufferedReader(new InputStreamReader( new BufferedInputStream(new FileInputStream(new File(GFF))),"utf-8"),5*1024*1024);// reading file using 5M buffer  	
	
		
	BufferedWriter bwAnoMap = new BufferedWriter(new OutputStreamWriter( new BufferedOutputStream(new FileOutputStream(new File(AnoMap))),"utf-8"),5*1024*1024);// Writing file using 5M buffer 
	
	int Sn=498; //Maximum Number of scaffolds
	int [] Ln = new int [Sn];//Maximum Number of nucleotides in a scaffold
	int Gn=2000000; //Maximum Number of genes
	String [][] RefGenome=new String [Sn][];
	int [][] AnoGene=new int [Sn][];
	int [][] AnoType=new int [Sn][];
	String [][] AnoStrand = new String [Sn][] ; //Gene strand: +1: sense; -1: antisense
	int [] Gene = new int[Gn] ; //Gene name
	int [] Type = new int[Gn] ; //Gene name
	int [] Sca = new int[Gn] ; //Scaffold
	int [] Start = new int[Gn] ; //Gene start position
	int [] End = new int[Gn] ; //Gene end position
	String [] Strand = new String[Gn] ; //Gene strand: +1: sense; -1: antisense
		
		
	BufferedReader brHead = new BufferedReader(new InputStreamReader( new BufferedInputStream(new FileInputStream(new File(HeadFile))),"utf-8"),16*1024);// reading file using 16k buffer 
	
	System.out.print("Reading the head file:\n"+HeadFile+"\n");
	System.out.print("Sca\t");
	System.out.print("Length\n");
	String lineHead="";
	String [] a = new String [50];
	int i=0;
	int j=0;
	int n=0; 
	int S=0; //Scaffold Number
	int L=0; //Locus Position	

	while((lineHead=brHead.readLine()) != null)
	{
		a = lineHead.replace(" ", "\t").split("\t");
		
		//skip the first/last/empty line
		if (lineHead.indexOf("scaffold_")>=0) 
		{					
			S=Integer.parseInt(a[0].substring(9)); //Scaffold
			Ln[S]= Integer.parseInt(a[1])+1; //Length of the scaffold
			RefGenome[S]=new String [Ln[S]];
			AnoGene[S]=new int [Ln[S]];
			AnoType[S]=new int [Ln[S]];
			AnoStrand[S] = new String [Ln[S]] ; 
			System.out.print(S+"\t");
			System.out.print(Ln[S]+"\n");
		}
	}
	brHead.close();
 	
	System.out.print("Reading the GFF file\n");
	
	System.out.print("Gene\tType\tSca\tStart\tEnd\tStrand\n");
	
	int NumScaffold=0; //Scaffold
	String lineGFF;
	String[] b = new String[20];	
	j=0; //Number of Annotations
	
	while((lineGFF=brGFF.readLine()) != null)
	{
		j++; 
		b = lineGFF.replace(" ", "\t").split("\t");
		Sca[j]= Integer.parseInt(b[0].substring(9)); //Scaffold
		
		int GeneType=0;
		if (b[2].equals("CDS")) {GeneType=1;} 
		if (b[2].equals("exon")) {GeneType=2;}
		if (b[2].equals("mRNA")) {GeneType=3;}
		if (b[2].equals("three_prime_UTR")) {GeneType=4;}
		if (b[2].equals("five_prime_UTR")) {GeneType=5;}
		
		Type[j]=GeneType;
		
		Start[j]= Integer.parseInt(b[3]); //Position
		End[j]= Integer.parseInt(b[4]); //Position
		Strand[j]=b[6];	
		int s=b[8].indexOf("gene")+4;
		int e=b[8].indexOf("-");
		if (s>1 && e>5)
		{
			Gene[j]=Integer.parseInt(b[8].substring(s,e));
			System.out.print("Gene"+Gene[j]+"\tScaffold_"+Sca[j]+"\t"+Start[j]+"-"+End[j]+"\t"+Strand[j]+"\t"+Type[j]+"-"+b[2]+"\n");
		}
	}
	
	Gn=j;
	
	System.out.print("Reading the Reference Genome:\n");

	String LineRef;
	S=0;
	int WrongChar=0;
	int UnkownBases_n=0;
	int UnkownBases_N=0;
	L=0;
	while((LineRef=brRefGenome.readLine()) != null)
	{
		int LL=0;
		
		LineRef=LineRef.replace(" ","").replace("\n","").replace("\r","").replace("\t","");

		int Ids=LineRef.indexOf(">scaffold_");
		
		if (Ids>=0) //the first line of a scaffold
		{
			S = Integer.parseInt(LineRef.substring(Ids+10)); //Scaffold
			L=0;
			//bwRefGenomeCopy.write("\n");			
			System.out.print("\n");
			//bwRefGenomeCopy.write(">scaffold_"+S+"\n");
			System.out.print("Reading Scaffold-"+S+"\n");
			
			if (S>NumScaffold) {NumScaffold = S;}
		}
		else
		{
			LL=LineRef.length();
			for (int k=0;k<LL;k++)
			{
				String Nucleotide=LineRef.substring(k,k+1);
				if(Nucleotide.equals("N"))
				{
					UnkownBases_N++;
				}
				if(Nucleotide.equals("n"))
				{
					UnkownBases_n++;
					//Nucleotide="N";
				}
				if("actgnACTGN".indexOf(Nucleotide)>=0)
				{
					L++;
					RefGenome[S][L]=Nucleotide;
					//bwRefGenomeCopy.write(RefGenome[S][L]);
					//System.out.print(RefGenome[S][L]);
				}
				else
				{
					WrongChar++;
				}
			}
		}
	}
	
	System.out.print("Number of Unkown Bases (n):"+UnkownBases_n+"\n");
	System.out.print("Number of Unkown Bases (N):"+UnkownBases_N+"\n");
	System.out.print("Number of Wrong Characters:"+WrongChar+"\n");
		
	//bwRefGenomeCopy.close();
		
		for (j=1;j<Gn;j++)
		{
			int s=Sca[j];
			for(int l=Start[j];l<=End[j];l++)
			{
				AnoGene[s][l]=Gene[j];
				AnoType[s][l]=Type[j];
				AnoStrand[s][l]=Strand[j];
			}
		}
			
	System.out.print("Write to file:\n"+AnoMap+"\n");
	
	String HeadLine="Sca\tPos\tRef\tGene\tType\tStrand\n";
	bwAnoMap.write(HeadLine);
	
	int TotalSites=0;
	for (int s=1;s<Sn;s++)//Scaffold
	{		
	  for (int l=1;l<Ln[s];l++)//Position
	  {
		TotalSites++;
		String GeneStr=AnoGene[s][l]>0?String.valueOf(AnoGene[s][l]):"";
		String TypeStr="";
		if (AnoType[s][l]==0) {TypeStr="";} 
		if (AnoType[s][l]==1) {TypeStr="CDS";} 
		if (AnoType[s][l]==2) {TypeStr="exon";}
		if (AnoType[s][l]==3) {TypeStr="intron";}
		if (AnoType[s][l]==4) {TypeStr="UTR3";}
		if (AnoType[s][l]==5) {TypeStr="UTR5";}
		if (AnoStrand[s][l]==null) {AnoStrand[s][l]="";} 
		String OutputStr=s+"\t"+l+"\t"+RefGenome[s][l]+"\t"+GeneStr+"\t"+TypeStr+"\t"+AnoStrand[s][l]+"\n";
		bwAnoMap.write(OutputStr);
	  }
	}	  
	bwAnoMap.close();
	System.out.print("Total Sites:"+TotalSites+"\n");

}
catch(Exception e)
{
	System.out.println("\nThis program is to CombineRefGFF\n"); 
	System.out.println("\nUsage: Java -cp ./ CombineRefGFF\n"); 
	System.out.println(e);
	e.printStackTrace();
}
}

}
/*
======================================================================
CombineRefGFF

Usage: CombineRefGFF using Gene Frequency Data and gene annotation file:
 
Java -cp ./ CombineRefGFF <Map-GeneFreqData> <GFF> > output
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
