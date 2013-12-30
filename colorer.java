
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.struc.MDocument;
import chemaxon.struc.Molecule;
import java.awt.Color;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import com.google.gson.Gson;

class Pairs{
	private int a;
	private int b;
	private String color;
	public void setA(int a){this.a=a;}
	public void setB(int b){this.b=b;}
	public void setColor(String color){this.color=color;}
	public String getColor(){return color;}
	public int getA(){return a;}
	public int getB(){return b;}
	public Color getRGB() {
		return Color.decode("0x"+color);//a color from the hex key;
	}
}

class Data {
	//well since we are bothering passing a json, maybe is worth to create a decent data structure once for all
    private String smiles;
    private List<Pairs> pairs;
    private String[] hmap;
    public String getSmiles() { return smiles; }
    public List<Pairs> getPairs() { return pairs; }
    public String[] getHmap(){
    	if (hmap ==null)
    		genHmap();
    	return hmap;}
    

    public void setSmiles(String smiles) { this.smiles = smiles; }
    public void setPairs(List<Pairs> pairs) { this.pairs = pairs; }
    public void genHmap(){
    	List<String> hm = new ArrayList<String>();
    	for (int i=0;i<pairs.size(); i++){
    		//System.out.println(i);
    		hm.add(i,pairs.get(i).getColor());
    		}
    	setHmap(hm);
    }
    public void setHmap(List<String> hmap) {this.hmap =  hmap.toArray(new String[hmap.size()]);} //WTF!!!
    public String toString() {
        return String.format("smiles: %s, pairs %s", smiles, pairs);
    }
}

public class colorer {
	
	/**
	 * @param args
	 */
	
    private static Color shrinkColors(Color c){
			//marvin only allows 64 different colors on a molecule. so, let's flatten stuff!
			//the idea is to have a sort of histogram: 256/64 = 4 --> [0-4] = 0, [5-8]=1 etc etc...
			int oldcolor = c.getGreen();
			int delta = c.getGreen()%4;
			String newcolor = Integer.toHexString(oldcolor - delta);
			c= Color.decode("0x00"+newcolor+"00");
			return c;
		}
    	
	
	public static void getPic(String id, String width,String height, String val){
		//generate picture with marvin and returns filename

		MDocument mdoc = drawChemicalColors(val);
        String filename= id+".svg";
        String params = "svg:w" + width + "h"+ height +"setcolors,amap,ez,chiral_all";
        try {
			byte[] d3 = MolExporter.exportToBinFormat(mdoc, params);
			FileOutputStream fos = new FileOutputStream(id);
			fos.write(d3);
			fos.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	
    public static MDocument drawChemicalColors(String value) {
		MDocument doc = getColorPic(value);
		return doc;
	}

	/**
	 * @param smi
	 * @param value
	 * @return mdocument
	 * Draw a colored image based on the reactivity of the functional group.
	 * 
	 * Since Rdkit and Chemaxon index the molecules in a different way, and the smiles are also needed to
	 * all the matching functions inside the structure, we execute the coloring is passing from the server
	 * a string like:
	 * C[CH2:1][O:1][CH2:1]N_000000#006fff
	 * where the "_" separates the mapped smiles from the colormap.
	 * then, is as simply as reading the indexed smiles, and assigning the color to the atoms with matching mapping
	 */
	private static MDocument getColorPic(String value) {
		
		Molecule mol;
		Data data = new Gson().fromJson(value, Data.class);
		String smi=data.getSmiles();
		try {
			mol = MolImporter.importMol(smi);
		} catch (MolFormatException e) {
			// if the generation fails return an empty image
			return  null;
		}
		MDocument mdoc = new MDocument(mol);
		
		
//		Data data=data1[0];
		String[] h_map = data.getHmap();//value.split("#");
		boolean lotsofcolors=false;
		//First we set and count the different colors before preparing the palette,
		Set<Color> colorset=  new HashSet<Color>();
		ArrayList<Color> colorarray=  new ArrayList<Color>();
		
		if (h_map.length>=64)
			lotsofcolors=true;
		//initialize the set and the array with the black;
		Color c = Color.black;
		colorset.add(c);
		colorarray.add(c);
		
				
		for (String col: h_map){
			c = Color.decode("0x"+col);//a color from the hex key
			
			if (lotsofcolors){
				
				c= shrinkColors(c);
			}
			colorset.add(c);
			colorarray.add(c);
		}
		
		//since sets are really terrible in java, we change them to arrays.
		Color[] colorlist = colorset.toArray(new Color[0]);
		Color[] colorarray2 = colorarray.toArray(new Color[0]);
		//we create the palette;
		//colorlist
		for (int i=0; i< colorset.size(); i++){
			//int j =Arrays.asList(colorlist).indexOf(colorarray2[i]);
			mdoc.setAtomSetColorMode(i, MDocument.SETCOLOR_SPECIFIED);
			mdoc.setBondSetColorMode(i,MDocument.SETCOLOR_SPECIFIED);
			mdoc.setAtomSetRGB(i, colorlist[0].getRGB()); //black atoms
			mdoc.setBondSetRGB(i, colorlist[i].getRGB()); //colored bonds
		}
		
		for (int i=0; i< mol.getBondCount(); i++){
		//now, we loop trough the molecule, and we assign each bond
		//to a set, based on it's mapping if it had a key different than 0;
			//bond  between number 2 and 19 has color #3? then assign it to group 3.
		  int start_idx = mol.getBond(i).getAtom1().getAtomMap();
		  int end_idx=mol.getBond(i).getAtom2().getAtomMap();
		  for (int j=0; j<data.getPairs().size(); j++){
			Pairs pair=data.getPairs().get(j);

			int a = pair.getA();
			int b= pair.getB();
			boolean samebond = rightPair(a,b,start_idx,end_idx);
			if (samebond){
				int group=0;
				if (lotsofcolors){
					group = Arrays.asList(colorlist).indexOf(shrinkColors(pair.getRGB()));					
				}
				else{
					group = Arrays.asList(colorlist).indexOf(pair.getRGB());
				}
				mol.getBond(i).setSetSeq(group);
				System.out.println("a= "+a+" b= "+b +"color= "+pair.getColor());
			}

			
		  }
		}
		//if the molecule has less than 4 atoms add Hydrogens.
		if (mol.getAtomCount()<4)  	mol.addExplicitHydrogens(1);
        return mdoc;
		
	}
	

	
	private static boolean rightPair(int a, int b, int start_idx, int end_idx) {
		// TODO Terrible.. but..
		if ((a == start_idx || a== end_idx) && (b == start_idx || b== end_idx))
			return true;
		return false;
	}


	public static void main(String[] args) {
		if (args.length!=4){
		System.out.println("usage: colorer filename width height json\n json should be " +
				"indexed : EG:C[CH2:1][O:1][CH2:1]N " +
				"\n colorset follows the format: 000000#00ff00");}	
//		String json="{'smiles': '[OH:1][CH2:2][c:3]1[cH:8][cH:7][n:6][cH:5][cH:4]1', 'pairs': [{'a': 3, 'color': '000000', 'b': 4}, {'a': 1, 'color': 'FF0000', 'b': 2}, {'a': 6, 'color': 'FF0000', 'b': 7}, {'a': 8, 'color': 'FF0000', 'b': 7}, {'a': 2, 'color': '00FF00', 'b': 3}, {'a': 8, 'color': '000000', 'b': 3}, {'a': 5, 'color': 'FF0000', 'b': 6}, {'a': 4, 'color': '000000', 'b': 5}]}";
		
//			String json="{"
//				+ "'smiles' : '[C:1][C:2][C:3][C:4]=[C:5]',"
//				+ "'pairs' : [{'a':1, b: 2, color:'000000'},{'a':2, b: 3, color:'00ff00'},{'a':3, b: 4, color:'002200'},{'a':4, b: 5, color:'006600'}]"//,{'a':[2,3,'00ff00']}]"	
//				+ "}";//"'pairs' : [{'a':1,'b':2,'color':'000000']},{'a':2,'b':3,'color':'000000']}]"	
					//String[] m={"[C:1]=[C:2][C:3]",""};//0_0_000000#1_2_005500#2_3_00ff00"};
					//"[Na+:235].[Na+:234].[CH2:188]=[CH:187][CH2:186][CH:184]([CH3:185])[CH:182]([CH3:183])[CH2:181][CH:180]([OH:189])[CH:179]([OH:190])[CH:177]1[CH2:176][CH2:175][C:173]2([CH3:191])[O:174][C:168]3([CH3:192])[CH2:167][C:165]4([CH3:193])[O:166][CH:160]5[CH:159]=[CH:158][CH2:157][CH:155]6[O:156][CH:150]7[CH2:149][CH:147]8[O:148][CH:141]9[CH2:140][C:138]%10([CH3:197])[O:139][C:134]([CH3:199])([CH:132]%11[CH2:131][CH2:130][CH:128]%12[O:129][C:123]%13([CH3:200])[CH2:122][CH:120]%14[O:121][C:114]%15([CH3:202])[CH2:113][CH:111]%16[O:112][C:106]%17([CH3:205])[CH2:105][CH:104]([OH:206])[CH:103]([CH:101]%18[O:102][CH:96]%19[CH:95]([OH:208])[CH:94]([OH:209])[CH:93]([CH2:92][CH:91]([OH:210])[CH2:90][CH:89]([OH:211])[CH:87]%20[O:88][CH:82]%21[CH:81]([OH:213])[CH:80]([OH:214])[CH:79]([CH:77]%22[O:78][CH:72]%23[CH:71]([OH:217])[CH:69]%24[O:70][CH:64]%25[CH2:63][CH:61]%26[O:62][CH:57]([CH2:56][CH:55]([OH:224])[CH:54]([OH:225])[CH:52]%27[O:53][CH:47]%28[CH:48]([CH2:50][CH:51]%27[OH:226])[O:49][C:43]%27([CH3:228])[CH2:42][CH:40]%29[O:41][C:35]%30([CH3:229])[CH2:34][CH:33]([OH:230])[CH:31]%31[O:32][CH:27]([CH:25]([CH3:26])[CH:24]([OH:233])[CH:1]([CH3:0])[CH2:2][CH2:3][CH:4]([O:19][S:20](=[O:22])(=[O:21])[O-:23])[CH:5]([OH:18])[CH:6]([CH3:7])[CH2:8][CH:9]([OH:17])[C:10](=[CH2:11])[C:12]([CH3:16])=[CH:13][CH2:14][OH:15])[CH:28]([OH:232])[CH:29]([OH:231])[CH:30]%31[O:37][CH:36]%30[CH2:38][CH:39]%29[O:45][CH:44]%27[CH:46]%28[OH:227])[CH:58]([O:219][S:220](=[O:222])(=[O:221])[O-:223])[CH:59]([OH:218])[CH:60]%26[O:66][CH:65]%25[CH2:67][CH:68]%24[O:74][CH:73]%23[CH:75]([OH:216])[CH:76]%22[OH:215])[O:84][CH:83]%21[CH2:85][CH:86]%20[OH:212])[O:98][CH:97]%19[CH2:99][CH:100]%18[OH:207])[O:108][CH:107]%17[CH:109]([OH:204])[C:110]%16([CH3:203])[O:116][CH:115]%15[CH2:117][CH2:118][C:119]%14([CH3:201])[O:125][CH:124]%13[CH2:126][CH:127]%12[O:133]%11)[CH:135]([OH:198])[CH2:136][CH:137]%10[O:143][C:142]9([CH3:196])[CH2:144][CH2:145][C:146]8([CH3:195])[O:152][C:151]7([CH3:194])[CH2:153][CH:154]6[O:162][CH:161]5[CH2:163][CH:164]4[O:170][CH:169]3[CH2:171][CH:172]2[O:178]1", 
					//"000200#000200#000900#000600#000000#000000#000000#000000#000200#000000#000000#000000#000000#000000#000000#000000#000000#000000#000000#000000#000000#000000#000000#000000#000000#001500#001500#007100#000300#000400#007a00#006c00#006c00#000500#00b300#007800#007c00#007a00#00d400#007600#007200#007200#00ec00#007800#007c00#007600#000300#006d00#007000#007000#007700#000400#006e00#006d00#000400#000100#004e00#005700#000400#000300#006100#005c00#005700#007900#005f00#006000#006000#008f00#006200#006200#005f00#000400#006200#006400#006200#000500#000500#006400#006200#005a00#000500#000500#005b00#005b00#005a00#005d00#000400#005100#005100#000100#001f00#000100#003700#004f00#000300#000300#005b00#005800#004f00#007200#000500#005f00#005b00#008600#000600#00a700#009200#008b00#008600#000500#009900#008400#008400#00ff00#009300#008b00#008b00#007300#007300#007200#008d00#008d00#00e100#007000#006e00#006e00#00c700#006d00#006500#006500#006700#006700#007e00#006d00#007400#000400#00b100#007e00#007e00#007d00#00e000#008200#007600#007600#006e00#006e00#006400#008000#008000#00ba00#006300#006d00#006c00#00b000#005d00#005000#005000#001400#000000#000000#003d00#005100#005100#00b800#005500#006200#003d00#00e700#004e00#004b00#004b00#00ae00#004a00#004200#004600#006800#006800#004700#004700#000100#000000#000b00#000100#000100#000100#000100#000000#000000#000000#000000#000100#004200#006b00#006200#007900#006400#008300#00ac00#000400#007400#00a600#007e00#00c100#00a100#000500#00af00#000600#000500#000300#000300#000100#000100#000400#000500#000500#000500#000500#000400#000300#000400#000000#000000#000000#000000#000100#000400#000400#000300#00c500#009200#000500#000400#000300#000000#000000#000000"};
//			getPic("mt.svg","800","800", json);
			getPic(args[0],args[1],args[2],args[3]);
		
	}
	
}
