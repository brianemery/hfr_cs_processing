function  [gridd,Lat,Lon]=sbgrid(opt)
% SBGRID.M - create grid for UCSB HFR Array
% [gridd,Lat,Lon]=sbgrid(opt)
%
% Makes the grid for computing total vectors in the format specified
% in rad2tot2.m, which is 2 columns, long and lat. All inputs are in km.
%
% Opt:      Outputs:
%  1        Original 2 and 3 site grid, used thru Oct 1998, Jan 2003 thru Mar 2003
%  2        5 site grid, used  from Nov 1998 thru Dec 2002
%  3        Original Grid, plus east channel coverage, used starting Apr 2003
%  4        Original Grid, plus east channel coverage, plus 5 site (Lon,Lat empty)
%  5        MEGA option, outputs a big square grid in Lat, Lon (gridd empty)
%  6        8+ site option, goes way north
%
% Gridspace' is HARDWIRED  to 2km.
% This was done because this mfile would need to be modified to truely be
% flexible in grid spacing
%
% USEFUL EXAMPLE
% % Get the gridd without land points
% gridd = sbgrid(3);
% kk = sbgrid_mask(gridd);
% mm = setdiff(1:size(gridd,1),kk);
% gridd = gridd(mm,:);
% 

% Copyright (C) 1997-2010 Brian Emery
% 21Dec97 - Version 1.0
% 09Feb99 - expanded for 5 sites
% 08Feb00 - Fixed redundant grid points, added output of
% 				indecies to keep.
% 13Sep00 - Re fixed redundant gridd points so they are not calculated in the first place
% 18May03 - Combined sbgrid.m, sbgrid2.m, and sbgrid3.m into one file, with selectablility
% 15Apr08 - added 5 site to east channel, expanded to incl SLO sites, 

% ------------------------------------------------------------
% DEFINITIONS
% -------------------------------------------------------------

%2 km grid spacing
gridspace=2; 

% get lat long of origin from codar_sites.m
loc = codar_sites({'central'});

% ------------------------------------------------------------
%  DEFINE GRIDS
% ------------------------------------------------------------

% Original Grid 
eastmin=-46;
eastmax=46;
northmin=-50;
northmax=14;
gridd1=makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc);

% Part of 5 site grid 
eastmin=-86;
eastmax=-48;
northmin=-50;
northmax=56;
gridd2=makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc);

% 3rd gridd (also part of 5 site grid)
eastmin=-46;
eastmax=-42;
northmin=16;
northmax=56;
gridd3=makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc);

% Eastern grid
eastmin=48;
eastmax=88;
northmin=-50;
northmax=14;
gridd4=makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc);
% Customize this one ... get rid of land points and near land points
gridd4=gridd4([1:21 27:44 48:66 70:230 232:250 253:271 274:291 295:312 316:333 ...
    337:353 358:372 379:392 400:411 421:431 442:451 463:470 484:489 505:509],:);

% Far west, and north grid
eastmin=-128;%-86;
eastmax=-88;%-48;
northmin=-50;%58;%-50;
northmax=100;%56;
gridd5=makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc);

% 6th grid (north - near San Luis)
eastmin=-86;%-46;
eastmax=-42;
northmin=58;
northmax=100;
gridd6=makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc);

% 7th grid in km 
eastmin=-128;
eastmax=168;
northmin=-140;%-70;
northmax=-52;
gridd7=makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc);

% 8th grid  
eastmin=90;
eastmax=168;
northmin=-50;
northmax=-20;
gridd8=makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc);

% ---------------------
% CHECK PLOT
% ---------------------

% sbchan_map
% soCalBight_map
% plot(gridd1(:,1),gridd1(:,2),'b.')
% plot(gridd2(:,1),gridd2(:,2),'y.')
% plot(gridd3(:,1),gridd3(:,2),'g.')
% plot(gridd4(:,1),gridd4(:,2),'go')
% plot(gridd5(:,1),gridd5(:,2),'m.')
% plot(gridd6(:,1),gridd6(:,2),'ro')
% plot(gridd7(:,1),gridd7(:,2),'c.')
% plot(gridd8(:,1),gridd8(:,2),'ko')
%  axis([-121.7 -118 33 35.5])
%  
% keyboard


% --------------------
%  BUILD GRIDS
% ---------------------

% Paste each part together. 
% This done in such a way that grid points in old total files with only 
% the 3 site grid will have indicies compatible with
% the other grids. (this might not work on the NEW 5+ site grid??)
%
% Opt:      Outputs:
%  1        Original 2 and 3 site grid, used thru Oct 1998, Jan 2003 thru Mar 2003
%  2        5 site grid, used  from Nov 1998 thru Dec 2002
%  3        Original Grid, plus east channel coverage, used starting Apr 2003
%  4        5 site PLUS east channel grid (MEGA option)

if opt==1
    % 3 site grid (include in all grids)
    gridd=[gridd1']';
%     Lat=[Lat1']';
%     Lon=[Lon1']';
Lat=[]; Lon=[];
    
elseif opt==2
    % 5 site grid appended to 3 site gridd
    gridd=[gridd1' gridd2' gridd3']';
    Lon=[]; Lat=[];%[Lon,Lat]=meshgrid([Lon2(:,1)' Lon1(:,1)'],[Lat2(1,:)]);
    
elseif opt==3
    % add east channel coverage
    gridd=[gridd1' gridd4']';
    % note here the the Lon, Lat are for the
    % square grid including points over land.
    %Lon=[Lon1' Lon4']'; Lat=[Lat1' Lat4']';
       % deal with these later
    Lon=[]; Lat=[];

elseif opt==4
    % append the 5 site grid to the 3 site + east channel grid
    gridd=[gridd1' gridd4' gridd2' gridd3']';
    % deal with these later
    Lon=[]; Lat=[];

elseif opt==5
    % The top secret option, which outputs a gridded Lat Lon of the whole
    % entire coverage area. gridd is empty.
    % xlats=[Lat1(:,1)' Lat2(:,1)' Lat4(:,1)'];
    xlons=[Lon1(:,1)' Lon2(:,1)' Lon4(:,1)'];
    % plot(xlons,xlats,'bs')
    ylats=[Lat2(1,:)];
    % ylons=[Lon2(1,:)];
    % plot(ylons,ylats,'bs')
    gridd=[];
    [Lon,Lat]=meshgrid(xlons,ylats);

elseif opt==6
    % .. and save as 2 column matrix, lons, then lats
    gridd=[gridd1' gridd4' gridd2' gridd3' gridd5' gridd6' gridd7' gridd8']';
    
    % deal with these at some time when needed
    Lat=[]; Lon=[];
    
    % Specific to option 6:
    maskIdx=defineGridMaskPoints;
    
    % get rid of gridd points in the maskIdx Index array
    kdx=setdiff(find(gridd(:,1)),maskIdx);
    gridd=gridd(kdx,:);
    
    %keyboard
    
end

% TESTING CODE:
% keyboard
% 
% plotgrid
% hold on
% sbchan_map
% soCalBight_map
% axis([-121.7 -118 33 35.5])
% 
% % test for repetitions
% b=unique(gridd,'rows');
% if length(b) ~= size(gridd,1)
%     disp('error?')
% end


end
% --------------------------------------------------------------------
function gridd = makeGrid(eastmin,eastmax,northmin,northmax,gridspace,loc)
% MAKE GRID
% uses the km2lonlat.m tool from Mike Cook for bw compatibility with
% Archvie HFR data. 

east=eastmin:gridspace:eastmax;
north=northmin:gridspace:northmax;

% use cook's mfile to convert km to lat long
[lon, lat] = km2lonlat_old(loc.central(1),loc.central(2),east,north);

% make a grid out of the lat long data
[Lat, Lon]=meshgrid(lat,lon);
[m,n]=size(Lat);

% reshape ...
gridd=[reshape(Lon,m*n,1) reshape(Lat,m*n,1)];
    
end
% --------------------------------------------------------------------
function maskIdx=defineGridMaskPoints
%
% Here's how I got this: 
% maskIdx=sbgrid_mask(gridd);
% But! inpolygon.m causes errors around the islands
% which can be seen with this:
% plotgrid
% hold on
% sbchan_map
% soCalBight_map
% axis([-121.7 -118 33 35.5])
% plot(gridd(maskIdx,1),gridd(maskIdx,2),'c.')


   maskIdx=[ 21     22     23     24     25     26     27     28     29     30     31     32     35     36    ...
        37     38     39     40     41     42     43     44     45     46     47     68     69     70  ...
        71     72     73     74     75     76     85     86     87     88     89     90     91     92 ...
        93     94    109    119    120    121    122    123    132    133    134    135    136    137  ...
        138    139    140    141    153    154    155    156    157    158    179    180    181    182  ...
        183    184    185    186    201    202    203    225    226    227    228    229    230    231   ...
        250    271    272    273   1126   1166   1167   1170   1171   1172   1173   1174   1211   1212 ...
        1213   1214   1215   1216   1217   1218   1219   1220   1221   1222   1233   1234   1255   1256 ...
        1257   1258   1259   1260   1261   1262   1263   1264   1265   1266   1267   1268   1269   1279 ...
        1280   1281   1282   1283   1284   1285   1286   1287   1288   1289   1296   1297   1298   1299  ...
        1300   1301   1302   1303   1304   1305   1306   1307   1308   1309   1310   1311   1312   1313   1314 ...
1315   1316   1326   1327   1328   1329   1330   1331   1332   1333   1334   1335   1336   1337   1338  ...
1339   1340   1341   1342   1343   1344   1345   1346   1347   1348   1349   1350   1351   1352   1353  ...
1354   1355   1356   1357   1358   1359   1360   1361   1362   1363   1372   1373   1374   1375   1376  ...
1377   1378   1379   1380   1381   1382   1383   1384   1385   1386   1387   1388   1389   1390   1391  ...
1392   1393   1394   1395   1396   1397   1398   1399   1400   1401   1402   1403   1404   1405   1406  ...
1407   1408   1409   1410   1419   1420   1421   1422   1423   1424   1425   1426   1427   1428   1429  ...
1430   1431   1432   1433   1434   1435   1436   1437   1438   1439   1440   1441   1442   1443   1444  ...
1445   1446   1447   1448   1449   1450   1451   1452   1453   1454   1455   1456   1457   1464   1465  ...
1466   1467   1468   1469   1470   1471   1472   1473   1474   1475   1476   1477   1478   1479   1480  ...
1481   1482   1483   1484   1485   1486   1487   1488   1489   1490   1491   1492   1493   1494   1495  ...
1496   1497   1498   1499   1500   1501   1502   1503   1504   1507   1508   1509   1510   1511   1512  ...
1513   1514   1515   1516   1517   1518   1519   1520   1521   1522   1523   1524   1525   1526   1527  ...
1528   1529   1530   1531   1532   1533   1534   1535   1536   1537   1538   1539   1540   1541   1542  ...
1543   1544   1545   1546   1547   1548   1549   1550   1551   1552   1553   1554   1590   1609   1629  ...
1955   1959   1960   3042   3043   3046   3049   3052   3070   3073   3076   3094   3095   3096   3097  ...
3099   3100   3102   3103   4721   4722   4744   4745   4768   4791   4814   4837   4860   4883   4906 ...
4929   4952   4973   4974   4975   4991   4992   4995   4996   4997   4998   5013   5014   5015   5016  ...
5017   5018   5019   5020   5021   5034   5035   5036   5037   5038   5039   5040   5041   5042   5043   ...
5044   5056   5057   5058   5059   5060   5061   5062   5063   5064   5065   5066   5067   5078   5079  ...
5080   5081   5082   5083   5084   5085   5086   5087   5088   5089   5090   5101   5102   5103   5104  ...
5105   5106   5107   5108   5109   5110   5111   5112   5113   5124   5125   5126   5127   5128   5129  ...
5130   5131   5132   5133   5134   5135   5136   5148   5149   5150   5151   5152   5153   5154   5155  ...
5156   5157   5158   5159   5172   5173   5174   5175   5176   5177   5178   5179   5180   5181   5182 ...
5194   5196   5197   5198   5199   5200   5201   5202   5203   5204   5205   5745   5746   5747   5748 ...
5749   5750   5894   5895   5896   5897   5898   5899   6043   6044   6045   6046   6047   6193   6194 ...
9973  10119  10120  10121  10122  10268  10269  10270  10271  10417  10418  10419  10420  10567  10568  ...
10569  10716  10717  10718  10865  10866  10867  11014  11015  11016  11162  11163  11164  11165  11311 ...
11312  11313  11314  11376  11377  11378  11379  11380  11381  11382  11383  11384  11385  11386  11460 ...
11461  11462  11463  11525  11526  11527  11528  11529  11530  11531  11532  11533  11534  11535  11608 ...
11609  11610  11611  11612  11674  11675  11676  11677  11678  11679  11680  11681  11682  11683  11684 ...
11685  11686  11757  11758  11759  11760  11761  11823  11824  11825  11826  11827  11828  11829  11830 ...
11831  11832  11833  11834  11837  11838  11839  11840  11841  11842  11843  11844  11845  11846  11847 ...
11848  11849  11850  11851  11852  11905  11906  11907  11908  11909  11910  11944  11945  11946  11947  ...
11948  11949  11950  11984  11985  11986  11987  11988  11989  11990  12008  12009  12023  12024  12025  ...
12026  12027  12028  12029  12030  12047  12048  12049  12050  12051  12052  12053  12054  12061  12062 ...
12063  12064  12065  12066  12067  12068  12069  12070  12082  12083  12084  12085  12086  12087  12088 ...
12089  12090  12091  12092  12093  12094  12095  12096  12097  12098  12099  12100  12101  12102  12103  ...
12104  12105  12106  12107  12108  12109  12110  12119  12120  12121  12122  12123  12124  12125  12126 ...
12127  12128  12129  12130  12131  12132  12133  12134  12135  12136  12137  12138  12139  12140  12141 ...
12142  12143  12144  12145  12146  12147  12148  12149  12150  12157  12158  12159  12160  12161  12162 ...
12163  12164  12165  12166  12167  12168  12169  12170  12171  12172  12173  12174  12175  12176  12177 ...
12178  12179  12180  12181  12182  12183  12184  12185  12186  12187  12188  12189  12190  12194  12195 ...
12196  12197  12198  12199  12200  12201  12202  12203  12204  12205  12206  12207  12208  12209  12210 ...
12211  12212  12213  12214  12215  12216  12217  12218  12219  12220  12221  12222  12223  12224  12225 ...
12226  12227  12228  12229  12230  12232  12233  12234  12235  12236  12237  12238  12239  12240  12241 ... 
12242  12243  12244  12245  12246  12247  12248  12249  12250  12251  12252  12253  12254  12255  12256 ...
12257  12258  12259  12260  12261  12262  12263  12264  12265  12266  12267  12268  12269  12270  12271 ...
12272  12273  12274  12275  12276  12277  12278  12279  12280  12281  12282  12283  12284  12285  12286 ...
12287  12288  12289  12290  12291  12292  12293  12294  12295  12296  12297  12298  12299  12300  12301 ...
12302  12303  12304  12305  12306  12307  12308  12309  12310  12311  12312  12313  12314  12315  12316 ...
12317  12318  12319  12320  12321  12322  12323  12324  12325  12326  12327  12328  12329  12330  12331 ...
12332  12333  12334  12335  12336  12337  12338  12339  12340  12341  12342  12343  12344  12345  12346 ...
12347  12348  12349  12350  12351  12352  12353  12354  12355  12356  12357  12358  12359  12360  12361 ...
12362  12363  12364  12365  12366  12367  12368  12369  12370  12371  12372  12373  12374  12375  12376 ...
12377  12378  12379  12380  12381  12382  12383  12384  12385  12386  12387  12388  12389  12390  12391 ...
12392  12393  12394  12395  12396  12397  12398  12399  12400  12401  12402  12403  12404  12405  12406 ...
12407  12408  12409  12410  12411  12412  12413  12414  12415  12416  12417  12418  12419  12420  12421 ...
12422  12423  12424  12425  12426  12427  12428  12429  12430  12431  12432  12433  12434  12435  12436 ...
12437  12438  12439  12440  12441  12442  12443  12444  12445  12446  12447  12448  12449  12450  12451 ...
12452  12453  12454  12455  12456  12457  12458  12459  12460  12461  12462  12463  12464  12465  12466 ...
12467  12468  12469  12470  12471  12472  12473  12474  12475  12476  12477  12478  12479  12480  12481 ...
12482  12483  12484  12485  12486  12487  12488  12489  12490  12491  12492  12493  12494  12495  12496  ...
12497  12498  12499  12500  12501  12502  12503  12504  12505  12506  12507  12508  12509  12510  12511 ...
12512  12513  12514  12515  12516  12517  12518  12519  12520  12521  12522  12523  12524  12525  12526 ...
12527  12528  12529  12530  12531  12532  12533  12534  12535  12536  12537  12538  12539  12540  12541 ...
12542  12543  12544  12545  12546  12547  12548  12549  12550];

% so... using this:
% nKeep=[]; while 1, n=gridnum(gridd), nKeep=[nKeep n]; end
% here are points near the islands to keep (ie, remove from the above
% list: 
keep=[21  11823  11824  11674  11675  11525  11526  11527  11376  11377  11378  11379  11385  11386 ...
    11535  11837     36     37  11837  11838  11839  11840  11849  11850  11851  11852   1554  ...
    1553   1552     47  11848  11534];

maskIdx=setdiff(maskIdx,keep);

end
