export ℎ, ħ, ℏ, ℯ, Φ₀, ϕ₀, ln2,
    sigma0, σ0, sigmax, σx, sigmay, σy, sigmaz, σz, sigmaplus, σ₊, sigmaminus, σ₋,
    ⊗

# Physical and mathematical constants
const ℎ = 6.626069934E-34 # Planck constant, in Joules*seconds, NIST 2017
const ħ = 1.054571783E-34 # Reduced Planck constant
const ℏ = ħ # Actual unicode symbol
const ℯ = 1.6021766208E-19 # Elementary charge, in Coulomb, 2014 CODATA
const Φ₀ = 2.067833831E-15 # Magnetic flux quantum (ℎ/2ℯ), in Weber, 2014 CODATA
const ϕ₀ = 3.291059757E-16 # Reduced flux quantum (ℏ/2ℯ), in Weber, 2014 CODATA
const ln2 = 0.693147180559945309 # Natural logarithm of 2
const sqrtln16 = 1.66510922231539551 # Square root of ln(16)

# Pauli matrices (spin-1/2 operators)
const sigma0 = qeye(2)
const sigmax = Operator(sparse([1,2],[2,1],[1.0,1.0],2,2),(2,),true)
const sigmay = Operator(sparse([1,2],[2,1],[-1.0im,1.0im],2,2),(2,),true)
const sigmaz = Operator(sparse([1,2],[1,2],[1.0,-1.0],2,2),(2,),true)
const σ0 = sigma0
const σx = sigmax
const σy = sigmay
const σz = sigmaz
const sigmaplus  = Operator(sparse([2],[1],[1.0],2,2),(2,),false) # σ₊*[1,0] = [0,1]
const sigmaminus = Operator(sparse([1],[2],[1.0],2,2),(2,),false)
const σ₊ = sigmaplus
const σ₋ = sigmaminus

# Function aliases
const ⊗ = kron
const tensor = kron
const fock = basis
const qzeros = qzero

# Table of square root of factorials √(n!) for n ∈ 0:127 (1KiB)
const sqrtfact_table = Float64[
    1.0000000000000000000, 1.0000000000000000000, 1.4142135623730950488,
    2.4494897427831780982, 4.8989794855663561964, 10.954451150103322269,
    26.832815729997476357, 70.992957397195392511, 200.79840636817813151,
    602.39521910453439454, 1904.9409439665052252, 6317.9743589223278835,
    21886.105181141755630, 78911.474450804681438, 295259.70128007648582,
    1.1435359058639129588E6, 4.5741436234556518350E6,
    1.8859677306253148084E7, 8.0014834285449844953E7,
    3.4877657663442939413E8, 1.5597762686284978865E9,
    7.1477928181858656891E9, 3.3526120082371710076E10,
    1.6078562354540587669E11, 7.8768547132293829008E11,
    3.9384273566146914504E12, 2.0082117944245961319E13,
    1.0434974580907397764E14, 5.5216695356722848749E14,
    2.9735100460129106405E15, 1.6286585271694955843E16,
    9.0679869067935485290E16, 5.1296280268036349125E17,
    2.9467469553410734788E18, 1.7182339742875652406E19,
    1.0165209277917570217E20, 6.0991255667505421303E20,
    3.7099532465014090857E21, 2.2869687743093501008E22,
    1.4282115417961529469E23, 9.0328029052332240864E23,
    5.7838159214452708158E24, 3.7483411234209726053E25,
    2.4579516484946125896E26, 1.6304206741784307882E27,
    1.0937194378152021970E28, 7.4179661362209580728E28,
    5.0855013667402369566E29, 3.5233386996620226535E30,
    2.4663370897634158575E31, 1.7439636808636061170E32,
    1.2454391808865586995E33, 8.9809896543167155890E33,
    6.5382591597917144388E34, 4.8046196242703894246E35,
    3.5632012788584194610E36, 2.6664556771205919519E37,
    2.0131298891248228833E38, 1.5331540468207617694E39,
    1.1776379687564843276E40, 9.1219444817107882595E40,
    7.1244663931920178495E41, 5.6098104478126475754E42,
    4.4526490041372451123E43, 3.5621192033097960898E44,
    2.8718723147247460219E45, 2.3331200978034608324E46,
    1.9097411059666879697E47, 1.5748128594969087944E48,
    1.3081378078327271991E49, 1.0944666130115569586E50,
    9.2221396029764281460E50, 7.8252449403763769254E51,
    6.6858922078602825325E52, 5.7514219472399922436E53,
    4.9808775141931966697E54, 4.3422283469044442401E55,
    3.8102899106011063425E56, 3.3651569321810681095E57,
    2.9910169058002623210E58, 2.6752468492881886262E59,
    2.4077221643593697636E60, 2.1802851503903891306E61,
    1.9863343046226278804E62, 1.8205054612841328324E63,
    1.6784231035053557904E64, 1.5565055535934567420E65,
    1.4518117296604018405E66, 1.3619201234191322394E67,
    1.2848328747704294744E68, 1.2188994890809338170E69,
    1.1627560052213890684E70, 1.1152763807523813626E71,
    1.0755335917960171153E72, 1.0427685057848376926E73,
    1.0163650175128549387E74, 9.9583027412855333880E74,
    9.8077907646157562109E75, 9.7092175013660332284E76,
    9.6605494379949297313E77, 9.6605494379949297313E78,
    9.7087320283538361861E79, 9.8053387065593642395E80,
    9.9513319291872585027E81, 1.0148407138632952789E83,
    1.0399022830248479171E84, 1.0706449288791818466E85,
    1.1074837259283487726E86, 1.1509308491181513962E87,
    1.2016070835354182374E88, 1.2602561412358770729E89,
    1.3277622343967479571E90, 1.4051714689748984594E91,
    1.4937177607097713179E92, 1.5948541417547209534E93,
    1.7102905289724946877E94, 1.8420392733196266468E95,
    1.9924701154117019517E96, 2.1643765498993678241E97,
    2.3610560690520684628E98, 2.5864073371085860453E99,
    2.8450480708194446498E100, 3.1424583053452915085E101,
    3.4851548555301426772E102, 3.8809042007129483894E103,
    4.3389828034793201975E104, 4.8704961173170509835E105,
    5.4887703709097355030E106
]
