# Reference and source map

**Executable version:** 1.3.4  
**Role:** supporting source map for the SSDP authority set and derived manual. External sources support specific imported facts, correlations, standards, or comparison data; they do not replace the package's D1-D4 authority.

## Scientific and numerical sources

**[R1] Janz et al. (1975), molten chloride density data.**  
G. J. Janz, R. P. T. Tomkins, C. B. Allen, J. R. Downey, Jr., G. L. Gardner, U. Krebs, and S. K. Singer, "Molten Salts: Volume 4, Part 2, Chlorides and Mixtures - Electrical Conductance, Density, Viscosity, and Surface Tension Data," *Journal of Physical and Chemical Reference Data* **4**(4), 871-1178 (1975).  
DOI: https://doi.org/10.1063/1.555527  
Used for the legacy LiCl, NaCl, KCl endpoint/binary density relations and their source-range interpretation.

**[R2] Van Artsdalen and Yaffe (1955), KCl-LiCl and KCl-NaCl measurements.**  
E. R. Van Artsdalen and I. S. Yaffe, "Electrical Conductance and Density of Molten Salt Systems: KCl-LiCl, KCl-NaCl and KCl-KI," *Journal of Physical Chemistry* **59**(2), 118-127 (1955).  
DOI: https://doi.org/10.1021/j150524a007  
Used as provenance for legacy KCl-LiCl and KCl-NaCl density behavior summarized by the Janz compilation.

**[R3] Fellner, Votava, and Chrenkova-Paucirova (1980), ternary LiCl-NaCl-KCl data.**  
P. Fellner, I. Votava, and M. Chrenkova-Paucirova, "Density of the Molten System LiCl-NaCl-KCl," *Chemical Papers / Chemicke Zvesti* **34**(3), 330-334 (1980).  
Article record: https://www.chempap.org/index.php/chempap/article/view/1898  
Full paper: https://www.chempap.org/file_access.php?file=343a330.pdf  
Used only as external ternary validation/comparison evidence. The package does **not** implement the Fellner ternary equation as its default Li/Na/K model.

**[R4] Kohler (1960), symmetric ternary extrapolation provenance.**  
F. Kohler, "Zur Berechnung der thermodynamischen Daten eines ternaren Systems aus den zugehorigen binaren Systemen," *Monatshefte fuer Chemie* **91**, 738-740 (1960).  
DOI: https://doi.org/10.1007/BF00899814  
Used as methodological provenance for the historical symmetric pair-weighting construction. The exact package formula is defined in `NUMERICAL_METHOD.md` and is the normative implementation target.

**[R5] Redlich and Kister (1948), binary excess-property representation provenance.**  
O. Redlich and A. T. Kister, "Algebraic Representation of Thermodynamic Properties and the Classification of Solutions," *Industrial & Engineering Chemistry* **40**(2), 345-348 (1948).  
DOI: https://doi.org/10.1021/ie50458a036  
Used as general provenance for the binary excess-property form. The package's LiCl-NaCl reconstruction uses a single equimolar anchor and the exact formula stated in `NUMERICAL_METHOD.md`; no higher-order Redlich-Kister fit is implied.

**[R12] Janz (1988), critically evaluated molten-salt correlations.**  
G. J. Janz, "Thermodynamic and Transport Properties for Molten Salts: Correlation Equations for Critically Evaluated Density, Surface Tension, Electrical Conductance, and Viscosity Data," *Journal of Physical and Chemical Reference Data* **17**, Supplement 2 (1988).  
NIST supplement index: https://www.nist.gov/srd/journal-physical-and-chemical-reference-data-monographs-or-supplements  
PDF: https://srd.nist.gov/JPCRD/jpcrdS2Vol17.pdf  
Used for the RbCl, CsCl, BeCl2, MgCl2, CaCl2, SrCl2, and BaCl2 pure-melt density correlations and tabulated temperature ranges.

## Structure and serialization sources

**[R6] AFLOW rock-salt prototype.**  
M. J. Mehl, D. Hicks, C. Toher, O. Levy, R. M. Hanson, G. Hart, and S. Curtarolo, "The AFLOW Library of Crystallographic Prototypes: Part 1," *Computational Materials Science* **136**, S1-S828 (2017).  
DOI: https://doi.org/10.1016/j.commatsci.2017.01.017  
Rock-salt prototype `AB_cF8_225_a_b-001`: https://aflow.org/p/AB_cF8_225_a_b-001/  
Used as crystallographic provenance for the conventional B1/rock-salt base template. The auxiliary chloride positions for divalent compositions are package-defined construction sites, not an asserted equilibrium crystal prototype.

**[R7] VASP POSCAR format.**  
VASP Wiki, "POSCAR."  
https://vasp.at/wiki/POSCAR  
Used for the POSCAR serialization contract.

**[R8] Crystallographic Information File (CIF) standard.**  
S. R. Hall, F. H. Allen, and I. D. Brown, "The Crystallographic Information File (CIF): a New Standard Archive File for Crystallography," *Acta Crystallographica Section A* **47**, 655-685 (1991).  
DOI: https://doi.org/10.1107/S010876739101067X  
Used for CIF-format provenance.

**[R9] LAMMPS atom style and data-file format.**  
LAMMPS documentation, `atom_style`: https://docs.lammps.org/atom_style.html  
LAMMPS documentation, `read_data`: https://docs.lammps.org/read_data.html  
Used for the `atomic` and `charge` data-file serialization contract.

## Units and formula-mass sources

**[R10] SI definitions and the Avogadro constant.**  
Bureau International des Poids et Mesures, *The International System of Units (SI Brochure)*, 9th edition (2019; continuously maintained, 2026 update).  
DOI: https://doi.org/10.59161/AUEZ1291  
The Avogadro constant is exactly `6.02214076e23 mol^-1` in the SI. Used for conversion from molar volume to finite-cell volume.

**[R11] NIST Chemistry WebBook, SRD 69.**  
National Institute of Standards and Technology, NIST Chemistry WebBook, SRD 69.  
DOI: https://doi.org/10.18434/T4D303  
https://webbook.nist.gov/chemistry/  
Used as a formula/molecular-weight reference. The package deliberately freezes the legacy Li/Na/K/Cl numerical masses used by prior versions; those frozen constants are compatibility parameters rather than a claim to be the latest recommended atomic-weight table.

## Source-to-authority mapping

- D1 scientific meaning and validity: `SCIENTIFIC_METHOD.md`, supported primarily by [R1], [R2], [R3], [R12].
- D2 numerical method: `NUMERICAL_METHOD.md`, with provenance from [R1]-[R5], [R10]-[R12].
- D3 architecture: `ARCHITECTURE.md`; crystallographic/output dependencies are supported by [R6]-[R9].
- D4 public software contract: `SPECIFICATION.md`; output-format requirements route to [R7]-[R9].
- `LiCl_NaCl_KCl_cell_builder_manual.md` is a derived technical guide and reproduces the full reference list for standalone use.
