import java.io.*;


/**
 * 
 */
public class TrivialCaller {
    
    private static final String SAM_SEPARATOR = "\t";
    private static final String VCF_SEPARATOR = "\t";
    
    private static final int MIN_SV_LENGTH = 40;
    
    /**
     * VCF fields
     */
    private static final int DEFAULT_QUAL = 60;
    private static final String DEFAULT_FILTER = "PASS";
    private static final String DEFAULT_GT = "0/1";
    private static int idGenerator = -1;
    
    /**
     * Chromosome buffer
     */
    private static String CHROMOSOMES_DIR;
    private static String currentChromosome;  // ID
    private static StringBuilder sb;  // Sequence
    
    
    /**
     * Remark: the program prints the calls to STDOUT.
     *
     * @param args 
     * 1=directory that contains one $Z.fa$ file for every chromosome $Z$. These
     * files can be created from a single reference file by doing e.g.: 
     *
     * cat ref.fasta | awk '{ if (substr($0,1,1)==">") { filename=(substr($1,2) ".fa") } print $0 >> filename; close(filename) }'
     *
     */
    public static void main(String[] args) throws IOException {
        int i, p, q;
        int flags, quality, pos, nAlignments;
        double rate;
        String str, chrom, cigar, seq;
        BufferedReader br;
        int[] tmpArray = new int[3];
        
        final String SAM_FILE = args[0];
        CHROMOSOMES_DIR=args[1];
        
        sb = new StringBuilder();
        br = new BufferedReader(new FileReader(SAM_FILE));
        str=br.readLine(); nAlignments=0;
        while (str!=null) {
            nAlignments++;
            p=str.indexOf(SAM_SEPARATOR); q=str.indexOf(SAM_SEPARATOR,p+1);
            flags=Integer.parseInt(str.substring(p+1,q));
            p=q; q=str.indexOf(SAM_SEPARATOR,p+1);
            chrom=str.substring(p+1,q);
            if (string2contig(chrom)==-1) {
                str=br.readLine();
                continue;
            }
            p=q+1; q=str.indexOf(SAM_SEPARATOR,p+1);
            pos=Integer.parseInt(str.substring(p,q));
            p=q+1; q=str.indexOf(SAM_SEPARATOR,p+1);
            quality=Integer.parseInt(str.substring(p,q));
            p=q+1; q=str.indexOf(SAM_SEPARATOR,p+1);
            cigar=str.substring(p,q);
            p=str.indexOf(SAM_SEPARATOR,q+1);
            p=str.indexOf(SAM_SEPARATOR,p+1);
            p=str.indexOf(SAM_SEPARATOR,p+1);
            q=str.indexOf(SAM_SEPARATOR,p+1);
            seq=str.substring(p+1,q);           
            call(flags,chrom,pos,quality,cigar,seq,tmpArray);
            if (nAlignments%10000==0) System.err.println("Processed "+nAlignments+" alignments");
            str=br.readLine();
        }
        br.close();
    }
    
    
    /**
     * Remark: every alignment is used, including low MAPQ, secondary, and 
     * supplementary, since the procedure tries to maximize recall without 
     * caring about precision.
     *
     * @param position 1-based.
     */
    private static final void call(int flags, String chrom, int pos, int quality, String cigar, String seq, int[] tmpArray) throws IOException {
        char c;
        int i, j, p;
        int length, first, last, startB, startCigar, endCigar;
        int referencePosition, seqPosition;
        
        cigar_prefixClip(cigar,tmpArray);
        startB=tmpArray[0];  // 0-based
        startCigar=tmpArray[1];  // 0-based
        cigar_length(cigar,startCigar,tmpArray);
        endCigar=tmpArray[0];  // 0-based
        referencePosition=pos-1;  // Last consumed pos (1-based).
        seqPosition=startB-1;  // Last consumed pos (0-based).
        p=startCigar;
        for (i=startCigar; i<=endCigar; i++) {
            c=cigar.charAt(i);
            if (c=='M' || c=='X' || c=='=') {
                length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
                referencePosition+=length; seqPosition+=length;
            }
            else if (c=='D' || c=='N') {
                length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
                if (length>=MIN_SV_LENGTH) printDel(chrom,referencePosition+1,referencePosition+length);
                referencePosition+=length;
            }
            else if (c=='I') {
                length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
                if (length<MIN_SV_LENGTH) seqPosition+=length;
                else {
                    if (!seq.equals("*")) printIns(chrom,referencePosition,seq.substring(seqPosition+1,seqPosition+length+1));
                    seqPosition+=length;
                }
            }
        }
    }
    
    
    /**
     * @param first,last the deleted positions are [first..last], 1-based.
     */
    private static final void printDel(String chrom, int first, int last) throws IOException {
        System.out.println(chrom+VCF_SEPARATOR+(first-1)+VCF_SEPARATOR+(++idGenerator)+VCF_SEPARATOR+getSubstring(chrom,first-2,last-1)+VCF_SEPARATOR+getChar(chrom,first-2)+VCF_SEPARATOR+DEFAULT_QUAL+VCF_SEPARATOR+DEFAULT_FILTER+VCF_SEPARATOR+"SVTYPE=DEL;SVLEN="+(last-first+1)+VCF_SEPARATOR+"GT"+VCF_SEPARATOR+DEFAULT_GT);
    }
    
    
    /**
     * @param pos the characters in $seq$ are inserted to the right of $pos$
     * (1-based).
     */
    private static final void printIns(String chrom, int pos, String seq) throws IOException {
        final char c = getChar(chrom,pos-1);
        System.out.println(chrom+VCF_SEPARATOR+pos+VCF_SEPARATOR+(++idGenerator)+VCF_SEPARATOR+c+VCF_SEPARATOR+(c+seq)+VCF_SEPARATOR+DEFAULT_QUAL+VCF_SEPARATOR+DEFAULT_FILTER+VCF_SEPARATOR+"SVTYPE=INS;SVLEN="+seq.length()+VCF_SEPARATOR+"GT"+VCF_SEPARATOR+DEFAULT_GT);
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
    
    
    /**
     * @param out output array: 
     * 0: length of the soft clip at the beginning of $cigar$; 
     * 1: the first position of $cigar$ from which to continue reading it after 
     *    the soft clip;
     * 2: length of the soft plus hard clip at the beginning of $cigar$.
     */
    private static final void cigar_prefixClip(String cigar, int[] out) {
        char c;
        int i, p;
        int length, clip, softClip;
        final int cigarLength = cigar.length();
        
        p=0; i=0; clip=0; softClip=0;
        while (i<cigarLength) {
            c=cigar.charAt(i);
            if (!Character.isDigit(c)) {
                if (c=='H') {
                    clip+=Integer.parseInt(cigar.substring(p,i));
                    p=i+1;
                }
                else if (c=='S') {
                    length=Integer.parseInt(cigar.substring(p,i));
                    softClip=length; clip+=length;
                    p=i+1;
                    break;
                }
                else break;
            }
            i++;
        }
        out[0]=softClip; out[1]=p; out[2]=clip;
    }
    
    
    /**
     * @param out output array: 
     * 0: length of the soft plus hard clip at the end of $cigar$;
     * 1: the last position of $cigar$ before the clip.
     */
    private static final void cigar_suffixClip(String cigar, int[] out) {
        boolean found;
        char c;
        int i, iPrime;
        int clip;
        final int cigarLength = cigar.length();
        
        // First hard clip, if any.
        clip=0;
        c=cigar.charAt(cigarLength-1);
        if (c=='H') {
            // Hard clip
            i=cigarLength-2;
            while (i>=0) {
                c=cigar.charAt(i);
                if (!Character.isDigit(c)) break;
                i--;
            }
            clip+=Integer.parseInt(cigar.substring(i+1,cigarLength));
            // Following soft clip, if any.
            c=cigar.charAt(i);
            if (c=='S') {
                iPrime=i; i--;
                while (i>=0) {
                    c=cigar.charAt(i);
                    if (!Character.isDigit(c)) break;
                    i--;
                }
                clip+=Integer.parseInt(cigar.substring(i+1,iPrime));
            }
        }
        else if (c=='S') {
            // Soft clip
            i=cigarLength-2;
            while (i>=0) {
                c=cigar.charAt(i);
                if (!Character.isDigit(c)) break;
                i--;
            }
            clip+=Integer.parseInt(cigar.substring(i+1,cigarLength-1));
        }
        else i=cigarLength-1;
        out[0]=clip; out[1]=i;
    }
    
    
    /**
     * @param out output array: 
     * 0: the last position before any soft and hard clip at the end of $cigar$,
     *    if any (call this position $to$);
     * 1: length of the alignment encoded in $cigar[from..to]$, projected on the
     *    reference; 
     * 2: length of the alignment encoded in $cigar[from..to]$, projected on the
     *    read.
     */
    private static final void cigar_length(String cigar, int from, int[] out) {
        char c;
        int i, p;
        int length, lengthA, lengthB;
        final int cigarLength = cigar.length();
        
        p=from; lengthA=0; lengthB=0;
        for (i=from; i<cigarLength; i++) {
            c=cigar.charAt(i);
            if (c=='M' || c=='X' || c=='=') {
                length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
                lengthA+=length; lengthB+=length;
            }
            else if (c=='I') {
                length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
                lengthB+=length;
            }
            else if (c=='D' || c=='N') {
                length=Integer.parseInt(cigar.substring(p,i)); p=i+1;
                lengthA+=length;
            }
            else if (c=='S') break;
        }
        out[0]=p-1; out[1]=lengthA; out[2]=lengthB;
    }
    
    
    /**
     * @return -1 if $str$ does not represent a standard contig.
     */
    private static final int string2contig(String str) {
        char c;
        int length, out;
        
        out=-1;
        length=str.length();
        if (length>3 && str.substring(0,3).equalsIgnoreCase("chr")) { str=str.substring(3); length-=3; }
        if (length==1) {
            c=str.charAt(0);
            if (c=='X' || c=='x') return 23;
            else if (c=='Y' || c=='y') return 24;
            else if (c=='M' || c=='m') return 25;
            else if (c=='*') return -1;
            else return Integer.parseInt(str);
        }
        else if (length==2) {
            try { out=Integer.parseInt(str); return out; }
            catch (Exception e) { return -1; }
        }
        else return -1;
    }
    
}
