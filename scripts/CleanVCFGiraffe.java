import java.util.Arrays;
import java.io.*;


/**
 * Given a sorted VCF, the program:
 * - Tries to fill in the SVLEN field of every call with supported type, using
 *   other fields created by short-read SV callers. This is necessary for DELLY.
 *   Remark: some SVLEN occurrences might not get added.
 * - Tries to transform every DUP to an INS. Remark: some DUPs might not get
 *   translated (e.g. those that are too long).
 * - Tries to set the REF field of a symbolic DEL to a sequence. Remark: some
 *   symbolic DELs might not get a REF sequence.
 * - Ensures that REF and ALT have the same first character, and that such 
 *   character is the same as in the ref.
 */
public class CleanVCFGiraffe {
    /**
     * REF or ALT sequences longer than this are not filled in.
     */
    private static int MAX_SV_LENGTH;
    
    private static int SHORT_READ_LENGTH = 200;
    
    /**
     * Chromosome buffer
     */
    private static String CHROMOSOMES_DIR;
    private static String currentChromosome;  // ID
    private static StringBuilder sb;  // Sequence
    
    
    /**
     * @param args 
     * 1=directory that contains one $Z.fa$ file for every chromosome $Z$. These
     * files can be created from a single reference file by doing e.g.: 
     *
     * cat ref.fasta | awk '{ if (substr($0,1,1)==">") { filename=(substr($1,2) ".fa") } print $0 >> filename; close(filename) }'
     * 
     * 3: 1=outputs only non-symbolic INS, and DEL, for which cleaning was 
     * successful. Sets all calls to PASS. Excludes non-standard chromosomes.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        CHROMOSOMES_DIR=args[1];
        MAX_SV_LENGTH=Integer.parseInt(args[2]);
        final boolean FOR_VG = Integer.parseInt(args[3])==1;
        final String OUTPUT_VCF = args[4];
        
        boolean success;
        char c;
        int i, p, q;
        int row, pos, end, length;
        int nCalls, nSkipped_wrongChr, nSkipped_wrongType_minusOne, nSkipped_wrongType_one, nSkipped_wrongType_two, nSkipped_chrBoundary, nSkipped_svlen, nSkipped_dupins, nSkipped_symbolicDel, nSkipped_symbolicIns;
        String str, tmpString;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        
        nCalls=0; nSkipped_wrongChr=0; nSkipped_wrongType_minusOne=0; nSkipped_wrongType_one=0; nSkipped_wrongType_two=0; nSkipped_chrBoundary=0; nSkipped_svlen=0; nSkipped_dupins=0; nSkipped_symbolicDel=0; nSkipped_symbolicIns=0;
        currentChromosome=""; sb = new StringBuilder();
        br = new BufferedReader(new FileReader(INPUT_VCF));
        bw = new BufferedWriter(new FileWriter(OUTPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.substring(0,6).equalsIgnoreCase("#CHROM")) bw.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
            if (str.charAt(0)==COMMENT) {
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            nCalls++;
            tokens=str.split("\t");
            if (FOR_VG && !isStandardChromosome(tokens[0])) {
			    str=br.readLine();
                nSkipped_wrongChr++;
			    continue;
            }
            row=svType2Row(getField(tokens[7],SVTYPE_STR));
            if (FOR_VG) {   
                if (row==-1) {
                    nSkipped_wrongType_minusOne++;
			        str=br.readLine();
			        continue;
                }
            }
            else if (row==-1) {
                bw.write(str); bw.newLine();
				str=br.readLine();
				continue;
            }
            pos=Integer.parseInt(tokens[1]);
            getChar(tokens[0],0);  // Just for loading the chr
            
            // Ensuring SVLEN
            success=true;
            tmpString=getField(tokens[7],SVLEN_STR);
            if (tmpString==null) {
                if (row==0 || row==1 || row==2) {
                    // DEL, INV, DUP.
                    tmpString=getField(tokens[7],END_STR);
                    if (tmpString!=null) {
                        length=Integer.parseInt(tmpString)-pos;
                        tokens[7]+=";SVLEN="+length;
                    }
                    else {
                        System.err.println("WARNING: Not adding SVLEN since END field cannot be found: "+str);
                        nSkipped_svlen++;
                        success=false;
                    }
                }
                else if (row==3) {
                    // INS
                    tmpString=getField(tokens[7],"INSLEN");
                    if (tmpString!=null) {
                        length=Integer.parseInt(tmpString);
                        tokens[7]+=";SVLEN="+length;
                    }
                    else {
                        tmpString=getField(tokens[7],"SVINSLEN");
                        if (tmpString!=null) {
                            length=Integer.parseInt(tmpString);
                            tokens[7]+=";SVLEN="+length;
                        }
                        else {
                            tmpString=getField(tokens[7],"SVINSSEQ");
                            if (tmpString!=null) {
                                length=tmpString.length();
                                tokens[7]+=";SVLEN="+length;
                            }
                            else {
                                if (!(tokens[4].length()>=4 && tokens[4].substring(0,4).equalsIgnoreCase("<INS"))) {
                                    length=tokens[4].length()-1;
                                    tokens[7]+=";SVLEN="+length;
                                }
                                else {
                                    System.err.println("WARNING: Not adding SVLEN since it cannot be inferred: "+str);
                                    nSkipped_svlen++;
                                    success=false;
                                }
                            }
                        }
                    }
                }
            }
            
            // Ensuring END
            tmpString=getField(tokens[7],END_STR);     
            if (tmpString==null) {
                if (row==0 || row==1 || row==2) {
                    // DEL, INV, DUP.
                    tmpString=getField(tokens[7],SVLEN_STR);
                    if (tmpString!=null) {
                        end=pos+Integer.parseInt(tmpString);
                        tokens[7]+=";END="+end;
                    }
                    else if (tokens[3].length()>1 && tokens[3].charAt(0)!='<') {
                        end=pos+tokens[3].length()-1;
                        tokens[7]+=";END="+end;
                    }
                    else System.err.println("WARNING: Not adding END since SVLEN or REF field cannot be found: "+str);
                }
                else if (row==3) {
                    // INS
                    tokens[7]+=";END="+(pos+1);
                }
            }
            
            // DUP to INS
            if (row==2) {  // DUP
                tmpString=getField(tokens[7],END_STR);
                if (tmpString==null) {
                    System.err.println("WARNING: Not converting DUP to INS since the END field cannot be found: "+str);
                    nSkipped_dupins++;
                    success=false;
                }
                else {
                    end=Integer.parseInt(tmpString);
                    if (end<=pos) {
                        System.err.println("WARNING: Not converting DUP to an INS since END<=POS: "+str);
                        nSkipped_dupins++;
                        success=false;
                    }
                    else if (end-pos>MAX_SV_LENGTH) {
                        System.err.println("WARNING: Not converting DUP to an INS since it is too long: "+str);
                        nSkipped_dupins++;
                        success=false;
                    }
                    else {
                        tokens[3]=getSubstring(tokens[0],pos-1,pos-1);
                        tokens[4]=getSubstring(tokens[0],pos-1,end-1);
                        p=tokens[7].indexOf(SVTYPE_STR);
                        q=tokens[7].indexOf(";",p+1);
                        if (q<0) q=tokens[7].length();
                        tokens[7]=tokens[7].substring(0,p)+SVTYPE_STR+"=INS"+tokens[7].substring(q);
                        p=tokens[7].indexOf(END_STR);
                        q=tokens[7].indexOf(";",p+1);
                        if (q<0) q=tokens[7].length();
                        tokens[7]=tokens[7].substring(0,p)+END_STR+"="+pos+""+tokens[7].substring(q);
                    }
                }
            }
            
            // Ensuring sequence REF for every DEL (symbolic and non-symbolic).
            // Non-symbolic DELs need to be fixed because Sniffles may not use
            // the correct sequence.
            if (row==0) {
                tmpString=getField(tokens[7],END_STR);
                if (tmpString==null) {
                    System.err.println("WARNING: Not adding REF to symbolic DEL since END cannot be found: "+str);
                    nSkipped_symbolicDel++;
                    success=false;
                }
                else {
                    end=Integer.parseInt(tmpString);
                    if (end-pos>MAX_SV_LENGTH) {
                        System.err.println("WARNING: Not adding REF to symbolic DEL since it is too long: "+str);
                        nSkipped_symbolicDel++;
                        success=false;
                    }
                    else {
                        tokens[3]=getSubstring(tokens[0],pos-1,end-1);
                        tokens[4]=getSubstring(tokens[0],pos-1,pos-1);
                    }
                }
            }
            
            // Trying to fix symbolic INS (using GATK-SV's INFO fields).
            if (tokens[4].equalsIgnoreCase("<"+INS_STR+">") || tokens[4].equalsIgnoreCase("<"+INS_ME_STR+">") || tokens[4].equalsIgnoreCase("<"+INS_ME_ALU_STR+">") || tokens[4].equalsIgnoreCase("<"+INS_ME_LINE1_STR+">") || tokens[4].equalsIgnoreCase("<"+INS_ME_SVA_STR+">") || tokens[4].equalsIgnoreCase("<"+INS_UNK_STR+">") || tokens[4].equalsIgnoreCase("<"+INS_NOVEL_STR+">")) {
                tmpString=getField(tokens[7],INSSEQ_STR);
                if (tmpString==null) tmpString=getField(tokens[7],SVINSSEQ_STR);
                if (tmpString!=null) {
                    tokens[3]=getSubstring(tokens[0],pos-1,pos-1);
                    tokens[4]=Character.toLowerCase(tmpString.charAt(0))==Character.toLowerCase(tokens[3].charAt(0))?tmpString:tokens[3].charAt(0)+tmpString;
                }
                else {
                    if (success) {
                        nSkipped_symbolicIns++;
                        success=false;
                    }
                }
            }
            
            // Ensuring that non-symbolic REF and ALT have the same first char
            if ((row==0 || row==3) && tokens[4].charAt(0)!='<') {
                c=getChar(tokens[0],pos-1);
                tokens[3]=c+tokens[3].substring(1);
                tokens[4]=c+tokens[4].substring(1);
            }
            
            // Forcing PASS
            if (FOR_VG) tokens[6]="PASS";
            
            // Outputting
            if (!FOR_VG || success) {
                bw.write(tokens[0]);
                for (i=1; i<tokens.length; i++) bw.write("\t"+tokens[i]);
                bw.newLine();
            }
            str=br.readLine();
        }
        br.close(); bw.close();
        if (FOR_VG) {
            System.err.println("Skipped calls:");
            System.err.println("Wrong chromosome: "+nSkipped_wrongChr+" ("+((100.0*nSkipped_wrongChr)/nCalls)+"%)");
            System.err.println("Type -1: "+nSkipped_wrongType_minusOne+" ("+((100.0*nSkipped_wrongType_minusOne)/nCalls)+"%)");
            System.err.println("Type 1: "+nSkipped_wrongType_one+" ("+((100.0*nSkipped_wrongType_one)/nCalls)+"%)");
            System.err.println("Type 2: "+nSkipped_wrongType_two+" ("+((100.0*nSkipped_wrongType_two)/nCalls)+"%)");
            System.err.println("Too close to chr boundary: "+nSkipped_chrBoundary+" ("+((100.0*nSkipped_chrBoundary)/nCalls)+"%)");
            System.err.println("Wrong length: "+nSkipped_svlen+" ("+((100.0*nSkipped_svlen)/nCalls)+"%)");
            System.err.println("DUP->INS fails: "+nSkipped_dupins+" ("+((100.0*nSkipped_dupins)/nCalls)+"%)");
            System.err.println("Symbolic DEL: "+nSkipped_symbolicDel+" ("+((100.0*nSkipped_symbolicDel)/nCalls)+"%)");
            System.err.println("Symbolic INS: "+nSkipped_symbolicIns+" ("+((100.0*nSkipped_symbolicIns)/nCalls)+"%)");
        }
    }
    
    
    /**
     * Remark: the procedure ensures that $chr$ is loaded in global variable
     * $sb$.
     *
     * @return $chr[first..last]$, where $first$ and $last$ are zero-based;
     * @param first,last the procedure resets them if they are out of range.
     */
    private static final String getSubstring(String chr, int first, int last) throws IOException {
        String str;
        BufferedReader br;
        
        if (!chr.equalsIgnoreCase(currentChromosome)) {
            currentChromosome=chr;
            sb.delete(0,sb.length());
            br = new BufferedReader(new FileReader(CHROMOSOMES_DIR+"/"+chr+".fa"));
            str=br.readLine();  // Skipping FASTA header
            str=br.readLine();
            while (str!=null) {
                sb.append(str);
                str=br.readLine();
            }
            br.close();
        }
        if (first<0) first=0;
        if (last>sb.length()-1) last=sb.length()-1;
        return sb.substring(first,last+1);
    }
    
    
    /**
     * Like the above, but for a single character.
     */
    private static final char getChar(String chr, int first) throws IOException {
        String str;
        BufferedReader br;
        
        if (!chr.equalsIgnoreCase(currentChromosome)) {
            currentChromosome=chr;
            sb.delete(0,sb.length());
            br = new BufferedReader(new FileReader(CHROMOSOMES_DIR+"/"+chr+".fa"));
            str=br.readLine();  // Skipping FASTA header
            str=br.readLine();
            while (str!=null) {
                sb.append(str);
                str=br.readLine();
            }
            br.close();
        }
        if (first<0) first=0;
        if (first>sb.length()-1) first=sb.length()-1;
        return sb.charAt(first);
    }
    
    
	private static final int svType2Row(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return 0;
		else if (type.equalsIgnoreCase(INV_STR)) return 1;
        else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return 2;
        else if ( type.equalsIgnoreCase(INS_STR) ||
                  type.equalsIgnoreCase(INS_ME_STR) ||
                  type.equalsIgnoreCase(INS_NOVEL_STR)
                ) return 3;
		else return -1;
	}

    
	/**
	 * @return NULL if $field$ does not occur in $str$.
	 */
	private static final String getField(String str, String field) {
		final int FIELD_LENGTH = field.length()+1;
		int p = str.indexOf(field+"=");
		if (p<0) return null;
		if (field.equalsIgnoreCase(END_STR) || field.equalsIgnoreCase(INSSEQ_STR) || field.equalsIgnoreCase(SVINSSEQ_STR)) {
			while (p>0 && str.charAt(p-1)!=';') p=str.indexOf(field+"=",p+1);
			if (p<0) return null;
		}
		final int q = str.indexOf(SEPARATOR,p+FIELD_LENGTH);
		return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
	}    
    
    
    /**
     * @return TRUE iff $str$ represents a chromosome in 1..23,X,Y,M.
     */
    private static final boolean isStandardChromosome(String str) {
        boolean found;
        char c;
        int length;
        
        length=str.length();
        if (length>3 && str.substring(0,3).equalsIgnoreCase("chr")) { str=str.substring(3); length-=3; }
        if (length==1) {
            c=str.charAt(0);
            return Character.isDigit(c) || c=='X' || c=='x' || c=='Y' || c=='y' || c=='M' || c=='m';
        }
        else if (length==2) return Character.isDigit(str.charAt(0)) && Character.isDigit(str.charAt(1));
        else return false;
    }
    
    
	/**
	 * Basic constants
	 */
	public static final char COMMENT = '#';
	public static final String SEPARATOR = ";";
	public static final String PASS_STR = "PASS";
	public static final String PRECISE_STR = "PRECISE";
	public static final String IMPRECISE_STR = "IMPRECISE";
	public static final String END_STR = "END";
	public static final String SVTYPE_STR = "SVTYPE";
	public static final String SVLEN_STR = "SVLEN";
	public static final String CHR2_STR = "CHR2";
	public static final String CT_STR = "CT";
	public static final String CT_325_STR = "'3to5'";
	public static final String CT_523_STR = "'5to3'";
	public static final String CT_525_STR = "'5to5'";
	public static final String CT_323_STR = "'3to3'";
    
	/**
	 * SV types: labels used by callers.
	 */
	public static final String DEL_STR = "DEL";
	public static final String DEL_ME_STR = "DEL:ME";
	public static final String DEL_INV_STR = "DEL/INV";
	public static final String INS_STR = "INS";
	public static final String INS_ME_STR = "INS:ME";
    public static final String INS_ME_ALU_STR = "INS:ME:ALU";
    public static final String INS_ME_LINE1_STR = "INS:ME:LINE1";
    public static final String INS_ME_SVA_STR = "INS:ME:SVA";
    public static final String INS_UNK_STR = "INS:UNK";
    public static final String INS_NOVEL_STR = "INS:NOVEL";
    public static final String INSSEQ_STR = "INSSEQ";
    public static final String SVINSSEQ_STR = "SVINSSEQ";
	public static final String DUP_STR = "DUP";
	public static final String DUP_TANDEM_STR = "DUP:TANDEM";
	public static final String DUP_INT_STR = "DUP:INT";
	public static final String INV_STR = "INV";
	public static final String INV_DUP_STR = "INVDUP";
	public static final String CNV_STR = "CNV";
	public static final String BND_STR = "BND";
	public static final String TRA_STR = "TRA";
    
	/**
     *
	 */
	public static final byte TYPE_INSERTION = 1;
	public static final byte TYPE_DELETION = 2;
	public static final byte TYPE_DEL_INV = 3;
	public static final byte TYPE_INVERSION = 4;
	public static final byte TYPE_INV_DUP = 5;
	public static final byte TYPE_DUPLICATION = 6;
	public static final byte TYPE_CNV = 7;
	public static final byte TYPE_BREAKEND = 8;
	public static final byte TYPE_TRANSLOCATION = 9;
    
	/**
	 * Confidence intervals of positions.
	 *
	 * Remark: some callers report a standard deviation instead of a confidence
	 * interval. Sniffles reports additional interval information in its BEDPE
	 * output (an alternative to VCF), which our programs disregard for 
	 * simplicity. Some callers use CILEN to express a "confidence interval 
	 * around inserted/deleted material between breakends": we interpret CILEN
	 * exactly like CIEND, and we ignore it for insertions (since representing
	 * variable-length insertion strings complicates our code).
	 */
	public static final String CI_SEPARATOR = ",";
	public static final String CIPOS_STR = "CIPOS";
	public static final String CIEND_STR = "CIEND";
	public static final String STD_START1_STR = "STD_quant_start";
	public static final String STD_START2_STR = "STD_POS1";
	public static final String STD_END1_STR = "STD_quant_stop";
	public static final String STD_END2_STR = "STD_POS2";
	public static final String CILEN_STR = "CILEN";
	public static final int CIPOS_STR_LENGTH = CIPOS_STR.length();
	public static final int CIEND_STR_LENGTH = CIEND_STR.length();
	public static final int STD_START1_STR_LENGTH = STD_START1_STR.length();
	public static final int STD_START2_STR_LENGTH = STD_START2_STR.length();
	public static final int STD_END1_STR_LENGTH = STD_END1_STR.length();
	public static final int STD_END2_STR_LENGTH = STD_END2_STR.length();
	public static final int CILEN_STR_LENGTH = CILEN_STR.length();
    
}