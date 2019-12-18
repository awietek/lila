#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

err1 = np.array([52.026478562084201,24.429104990496157,13.784432903537898,8.6315229812780743,5.7969752975198148,3.9378916101394097,2.5583392271234828,1.7504987006372019,1.2372972699221236,0.94344927457903083,0.76867476461616491,0.66180463281435919,0.59811242725523783,0.50319245555345304,0.31493476649566787,0.19041978458508879,0.10995848379476314,0.063282768871992801,0.036179603573856411,0.020779306666106834,0.012394168423547569,0.0076065261031246223,0.0047175715519998107,0.0028099646757411278,0.0016583950596711361,0.00095883532178930864,0.00055488223324573482,0.00032057749010760972,0.00018752234041130578,0.00010633075582688889,5.9045958685999267e-05,3.3146251183779896e-05,1.9380569654003921e-05,1.1775436412619911e-05,7.426849307989869e-06,4.5404716431107772e-06,2.6330893447834569e-06,1.488673447624933e-06,8.0462731233410523e-07,4.1981778053923335e-07,2.1649962178571514e-07,1.1035913161094868e-07,6.0086684072757635e-08,3.3597814308450324e-08,1.9083415736531606e-08,1.0321848264993605e-08,5.9298699284227041e-09,3.4905198731394194e-09,2.0058195104866172e-09,1.1409539979467809e-09,6.2606630990558187e-10,3.5203129300498404e-10,2.014317601606308e-10,1.0868461686186492e-10,5.4839688345964532e-11,2.8343549729470396e-11,1.5894840998953441e-11,8.8107299234252423e-12,4.6469494918710552e-12,2.4726887204451486e-12,1.3002932064409833e-12,6.6080474425689317e-13,3.0553337637684308e-13,9.2370555648813024e-14,-3.5527136788005009e-14,-1.0658141036401503e-13,-1.4921397450962104e-13,-1.7763568394002505e-13,-1.9184653865522705e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-1.9184653865522705e-13,-1.9184653865522705e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-1.9184653865522705e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-2.1316282072803006e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.9184653865522705e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-1.8474111129762605e-13,-2.0605739337042905e-13,-1.9184653865522705e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.9184653865522705e-13,-2.0605739337042905e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.8474111129762605e-13,-1.9895196601282805e-13,-1.9184653865522705e-13,-1.9184653865522705e-13,-2.0605739337042905e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.9184653865522705e-13,-2.0605739337042905e-13,-2.0605739337042905e-13,-2.1316282072803006e-13,-2.0605739337042905e-13,-2.1316282072803006e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-1.9184653865522705e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-1.9184653865522705e-13,-2.0605739337042905e-13,-2.0605739337042905e-13,-1.9184653865522705e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-2.1316282072803006e-13,-2.0605739337042905e-13])


err2 = np.array([50.00826884767757,24.338117326076226,14.551360949282277,9.5669587371398634,6.0612833525930085,3.7816249085893361,2.7282125494961136,2.1134086924119657,1.6725880728904414,1.3053983014246953,1.0025731599546859,0.65816862535275789,0.29743494048238261,0.059545790744138571,0.033866896562855686,0.020414584809756775,0.01259202095904044,0.008291940005825893,0.0055249298871231645,0.0036989843285795132,0.002554576243788631,0.0017585845699485958,0.0012227054569393658,0.00081762953379893588,0.00052184767321250547,0.00032452583543118863,0.00019872972025325453,0.00011830960720260464,6.9178909662070964e-05,4.1977494198874865e-05,2.7480256889589327e-05,1.8739874064976902e-05,1.2583594504178564e-05,8.0650112934677054e-06,4.9884850668036052e-06,3.0449694037315567e-06,1.8934053613861579e-06,1.182833223367652e-06,7.4775064007326364e-07,4.9793846557122379e-07,3.2449849385329799e-07,1.9558688535425972e-07,1.14164961928509e-07,6.4388117948510626e-08,3.574314177967608e-08,1.9448435750746285e-08,1.0937441174974083e-08,6.3823151208453055e-09,3.8094540855126979e-09,2.3056685449773795e-09,1.4322978358904948e-09,9.071143836081319e-10,5.8879123798760702e-10,3.7837111221961095e-10,2.3743496058159508e-10,1.468407617721823e-10,8.7197804532479495e-11,5.1848303428414511e-11,2.9416469260468148e-11,1.7408297026122455e-11,1.0693668173189508e-11,6.8212102632969618e-12,4.4124703890702222e-12,2.8563817977556027e-12,1.8332002582610585e-12,1.1155520951433573e-12,6.2527760746888816e-13,3.1974423109204508e-13,1.0658141036401503e-13,-2.1316282072803006e-14,-9.9475983006414026e-14,-1.3500311979441904e-13,-1.7053025658242404e-13,-1.7763568394002505e-13,-1.9184653865522705e-13,-1.9895196601282805e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-2.1316282072803006e-13,-2.0605739337042905e-13,-2.1316282072803006e-13,-2.0605739337042905e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.0605739337042905e-13,-2.2737367544323206e-13,-2.0605739337042905e-13,-1.9895196601282805e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.2737367544323206e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.2737367544323206e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.2737367544323206e-13,-2.2737367544323206e-13,-2.2737367544323206e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-1.9895196601282805e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.0605739337042905e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.2026824808563106e-13,-2.2026824808563106e-13,-2.1316282072803006e-13,-2.1316282072803006e-13])

err3 = np.array([50.81523280847037,23.955477921348624,14.25489100612419,9.2180474401099772,6.2611846176882366,4.6918671508087542,3.5551491581890531,2.5728378920175246,1.7758942098600983,1.1822377026658089,0.77801056679388836,0.52617015792853294,0.38342838977461469,0.30480364679899452,0.25528435575728992,0.22402040617230057,0.20346992370608064,0.18782367179446879,0.17311266952857096,0.15631128775221015,0.13215942167247618,0.095107857071639046,0.064179654039740797,0.040539988784992431,0.026085415840803705,0.01743635027821,0.012430588467040593,0.0091567092832605113,0.0071428541407172474,0.0055960753002608499,0.0043404974419019027,0.0033798778733711288,0.0025944877983263837,0.0019581280260041467,0.001414673864950089,0.00099112557411018543,0.00068080944979698188,0.00047310761433294601,0.00032539247103358093,0.00022192450062163971,0.00015528540316722683,0.00010967163775177369,7.8365146016778908e-05,5.3416114838000794e-05,3.5612934539130947e-05,2.3908173666598032e-05,1.7207907134775269e-05,1.2755791054530619e-05,9.3625709354228093e-06,6.7379027015590509e-06,4.7156382834145916e-06,3.3020965162222637e-06,2.3761217633477827e-06,1.7632528326316788e-06,1.3279457391490723e-06,1.0198723074950067e-06,7.4204463373916951e-07,4.9214131081498635e-07,3.2191737631137585e-07,2.1478474110381285e-07,1.4281643956337575e-07,9.4101132219748251e-08,6.2816852164360171e-08,4.1799644634465949e-08,2.8265823459605599e-08,2.0099328423839324e-08,1.4778301249407377e-08,1.1192156534889364e-08,8.8731226810523367e-09,7.1223169584300194e-09,5.6458944186488225e-09,4.3036720853706356e-09,3.1494167274104257e-09,2.3214994371301145e-09,1.6985524098345195e-09,1.2344401056907373e-09,8.5874063415758428e-10,5.5718629710099776e-10,3.5018388189200778e-10,2.1758950197181548e-10,1.3996270809002453e-10,9.2519769623322645e-11,6.6116001562477322e-11,4.8842707656149287e-11,3.5633718198369024e-11,2.5508484213787597e-11,1.9234391857025912e-11,1.4814816040598089e-11,1.177369313154486e-11,9.4573238129669335e-12,7.595701845275471e-12,5.9330318435968366e-12,4.4408920985006262e-12,3.0624391911260318e-12,1.9539925233402755e-12,1.3216094885137863e-12,9.4502183856093325e-13,7.0343730840249918e-13,5.1159076974727213e-13,3.979039320256561e-13,3.2684965844964609e-13,2.8421709430404007e-13,2.4158453015843406e-13,2.2026824808563106e-13,1.9184653865522705e-13,1.8474111129762605e-13,1.6342482922482304e-13,1.5631940186722204e-13,1.5631940186722204e-13,1.4921397450962104e-13,1.4921397450962104e-13,1.4210854715202004e-13,1.4921397450962104e-13,1.3500311979441904e-13,1.3500311979441904e-13,1.3500311979441904e-13,1.2789769243681803e-13,1.2789769243681803e-13,1.3500311979441904e-13,1.3500311979441904e-13,1.3500311979441904e-13,1.3500311979441904e-13,1.3500311979441904e-13,1.2789769243681803e-13,1.4210854715202004e-13,1.4210854715202004e-13,1.3500311979441904e-13,1.4210854715202004e-13,1.2789769243681803e-13,1.3500311979441904e-13,1.2789769243681803e-13,1.3500311979441904e-13,1.4210854715202004e-13,1.3500311979441904e-13,1.4210854715202004e-13,1.3500311979441904e-13,1.3500311979441904e-13,1.2789769243681803e-13,1.4210854715202004e-13,1.4210854715202004e-13,1.3500311979441904e-13,1.4210854715202004e-13,1.4210854715202004e-13,1.3500311979441904e-13,1.3500311979441904e-13,1.3500311979441904e-13])

err4 = np.array([49.339219783070753,25.669030434201922,15.252368636630145,9.3556768905531058,6.4641394853568102,4.9274694677629469,3.9159826116091097,3.0330643104135788,2.3935624389707542,1.8608833480791205,1.4359625715538016,1.0754775951562081,0.82865262136618867,0.6655515644641099,0.55623260633878857,0.46308024050508578,0.37523079525635694,0.28545635235150257,0.20179832904821637,0.12361688466131682,0.061838947437486524,0.035802442816404323,0.024916105359835683,0.017651705873447554,0.01301675121910506,0.0097340695230059282,0.0074268754259207981,0.0058582433467222472,0.0047071965923635162,0.0038259559607638494,0.003157434453704866,0.0026306540342417861,0.00215707562808376,0.0017255082735729843,0.0013110660619730652,0.00096569683709191168,0.00070438263890792996,0.00051073298089931995,0.00037305204046589324,0.00028412602036809176,0.00021809370719694243,0.00016095276468774955,0.0001151016586149467,7.8709118440656312e-05,5.3837135581602524e-05,3.6629954983879998e-05,2.5808478461897266e-05,1.8898633868502657e-05,1.3669647970004917e-05,9.8467128495371981e-06,7.2435979703300291e-06,5.3156685382305113e-06,3.875481780823975e-06,2.7765893406694886e-06,1.9572201637174658e-06,1.385883514615216e-06,9.758431360751274e-07,6.737560269698406e-07,4.3785771453030975e-07,2.8272760488334825e-07,1.8672763246740942e-07,1.2560153095364512e-07,8.8828329580792342e-08,6.3367849634232698e-08,4.5212878774236742e-08,3.2959654561182106e-08,2.4909418527840899e-08,1.9195553591089265e-08,1.5136009778871085e-08,1.2170453089765942e-08,9.7363965778640704e-09,7.5583770353659929e-09,5.7275144627055852e-09,4.2723655724330456e-09,3.2275409012072487e-09,2.4956392508101999e-09,1.9448336274763278e-09,1.4981083040765952e-09,1.1255991694270051e-09,8.4201445815779152e-10,6.1821481267543277e-10,4.5717030161540606e-10,3.375575374775508e-10,2.4591173541921307e-10,1.7472956415076624e-10,1.2167333807155956e-10,8.3822726537619019e-11,6.0076388308516471e-11,4.5012882310402347e-11,3.5157654565409757e-11,2.8009594643663149e-11,2.1898927116126288e-11,1.6122214674396673e-11,1.0935252703347942e-11,6.8638428274425678e-12,4.4764192352886312e-12,3.1050717552716378e-12,2.1387336346379016e-12,1.4352963262354024e-12,9.9475983006414026e-13,7.3185901783290319e-13,5.8975047068088315e-13,4.7606363295926712e-13,3.765876499528531e-13,3.0553337637684308e-13,2.2737367544323206e-13,1.6342482922482304e-13,1.3500311979441904e-13,9.9475983006414026e-14,7.1054273576010019e-14,4.2632564145606011e-14,1.4210854715202004e-14,-1.4210854715202004e-14,-4.2632564145606011e-14,-5.6843418860808015e-14,-6.3948846218409017e-14,-6.3948846218409017e-14,-7.1054273576010019e-14,-7.815970093361102e-14,-8.5265128291212022e-14,-9.2370555648813024e-14,-8.5265128291212022e-14,-8.5265128291212022e-14,-7.815970093361102e-14,-8.5265128291212022e-14,-9.2370555648813024e-14,-8.5265128291212022e-14,-8.5265128291212022e-14,-9.2370555648813024e-14,-9.2370555648813024e-14,-1.0658141036401503e-13,-8.5265128291212022e-14,-9.2370555648813024e-14,-8.5265128291212022e-14,-9.2370555648813024e-14,-8.5265128291212022e-14,-8.5265128291212022e-14,-9.2370555648813024e-14,-8.5265128291212022e-14,-7.815970093361102e-14,-8.5265128291212022e-14,-8.5265128291212022e-14,-8.5265128291212022e-14,-8.5265128291212022e-14,-8.5265128291212022e-14,-8.5265128291212022e-14])

err5 = np.array([50.007063427494415,26.214974517372557,15.023794125223667,9.9363487079897013,6.9491805648964942,4.8811151131949089,3.6453664467781266,2.8756100643333014,2.2342541314694486,1.6897477690366003,1.2312318814911976,0.86434791709909575,0.62373522171294837,0.45394073318941963,0.34083856475265151,0.25665903504129517,0.19767748700514431,0.15678112012982837,0.13098191902589207,0.11117159706259372,0.09673707113326202,0.084888511757178264,0.075094455096760271,0.066337487704593912,0.059086452404301326,0.053300570930602476,0.048501834853482251,0.044850080071228149,0.042164390278912833,0.03999358273261322,0.038333605525764369,0.036957817127337478,0.035627723056464333,0.034325118399870291,0.032897467943413972,0.031431459586450217,0.029931461731372622,0.028221003494422803,0.026282915649154859,0.024149071689535617,0.021652011788660275,0.018614307494701166,0.015262523605642286,0.012000935572274329,0.0095890815651316075,0.0076728873500258032,0.0062162238088774302,0.005098458523860927,0.0041047981237340991,0.0033004687355600026,0.0027021772125905841,0.0022342469524190278,0.0018584518287809715,0.0015447607329690527,0.0012794153513695505,0.0010527891040581494,0.00085940802042472342,0.00068918081535684905,0.00053318079777397998,0.00041121777203301235,0.00033178585918136605,0.00027792318099528757,0.00024426139530220325,0.00021928355251077392,0.00020030416099103832,0.00018686254734490149,0.00017767361026699291,0.00017018821777980975,0.00016381561171385783,0.00015843577386931429,0.0001533944304128454,0.00014822487280952146,0.00014216346962570015,0.00013494036183914204,0.00012745269687286509,0.00012046733679227373,0.00011361162173528783,0.0001067935568954681,9.8952914598271491e-05,9.0223579327641801e-05,8.0559778446342989e-05,7.0731827761960631e-05,6.067060224523857e-05,5.0141556450000735e-05,4.011248650925836e-05,3.1712886553236785e-05,2.4217464194009608e-05,1.8921931953741478e-05,1.5002462504298819e-05,1.1959730599642171e-05,9.5568579752125515e-06,7.4590032355104086e-06,5.6682428066778812e-06,4.0931288509682417e-06,2.7948217606876824e-06,1.910703751661913e-06,1.3673505350197956e-06,9.5575705216788265e-07,6.565066925645624e-07,4.5401196047123449e-07,3.2714937958644441e-07,2.4892560190892254e-07,1.9652991767316053e-07,1.568701932797012e-07,1.2538535543171747e-07,1.0064799482734088e-07,8.1351885228286847e-08,6.8764734351134393e-08,5.9618734837840748e-08,5.2173916742503934e-08,4.491320737542992e-08,3.8376654742933169e-08,3.2434030572403572e-08,2.6910406347724347e-08,2.214196825889303e-08,1.8225605913357867e-08,1.4999407937921205e-08,1.2207237887196243e-08,9.8449888241702865e-09,8.1111437566505629e-09,6.9030932081659557e-09,6.1056866229591833e-09,5.4972559837551671e-09,4.9706940785654297e-09,4.4956252054362267e-09,4.067921111072792e-09,3.6689016269519925e-09,3.2435849561807117e-09,2.7456010798232455e-09,2.1784600789942488e-09,1.6603749486421293e-09,1.2885124078820809e-09,1.0254765925310494e-09,8.4118312315695221e-10,7.1688788239043788e-10,6.3444360876019346e-10,5.5563020850968314e-10,4.815916554434807e-10,4.2117420662179939e-10,3.7352521076172707e-10,3.3562486123628332e-10,3.0149749363772571e-10,2.6999913416148047e-10,2.3953816707944497e-10,2.134612486770493e-10,1.9334578382768086e-10])

plt.semilogy(err1, label="1")
plt.semilogy(err2, label="2")
plt.semilogy(err3, label="3")
plt.semilogy(err4, label="4")
plt.semilogy(err5, label="5")
plt.legend()
plt.show()