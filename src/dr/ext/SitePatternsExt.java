package dr.ext;

import alignment.Haplotypes;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.util.TaxonList;

public class SitePatternsExt extends SitePatterns{

	public SitePatternsExt(Alignment alignment, TaxonList taxa, int from,
			int to, int every, boolean strip) {
		super(alignment, taxa, from, to, every, strip);

	}
	public SitePatternsExt(Haplotypes haplotypes, TaxonList taxa, int from,
			int to, int every, boolean strip) {
//		Alignment alignment = haplotypes.getAlignment();
		this(haplotypes.getAlignment(), taxa, from, to, every, strip);
	}
	public void updateAlignment(Alignment alignment){
//		this.strip = strip;
//	    this.unique = unique;
	        
        setPatterns(alignment, from, to, every);
	}
	/**
	 * 
	 */
	private static final long serialVersionUID = 8376785905234747979L;

//	SitePatterns patterns2 = new SitePatterns(alignment, null, 0, -1, 1, true);
	
}
