# Pomáhá oktonionový asociátor v MPPT? Nezávislé numerické vyhodnocení

**Autor:** karagos01 · **Typ:** pracovní text / technický komentář · **Datum:** 25. 9. 2026 ·
**Anglická verze:** [PAPER.md](PAPER.md)
**Licence:** tento text CC BY 4.0 · přiložený kód MIT

**Komentář k:** M. Mazgal, *Intensional Maxwellian Formalism in Power Electronics: Causal Octonion
Control of Distributed Multi-Phase Converters*, v 1.0 (DOI 10.5281/zenodo.22877912) a v 1.1
(DOI 10.5281/zenodo.22914238), a k repozitáři `Causal-Octonion-MPPT.`

---

## Abstrakt

Nedávný preprint navrhuje řídit fotovoltaický MPPT převodník pomocí asociátoru neasociativní
oktonionové algebry, `[X, Y, Z] = (XY)Z − X(YZ)`, jako projekce „skryté kauzální paměti“ do střídy
PWM. Návrh neuvádí ani pravidla násobení té deformované algebry, ani fyzikální veličiny vstupující do
osmi složek oktonionu, ani převod asociátoru na střídu, a neobsahuje žádná měření. Tento komentář
chybějící definice doplňuje v nejvěrnější podobě, jakou jsme dokázali rekonstruovat, implementuje
kompletní oktonionový součin a vzniklé regulátory vyhodnocuje ve fotovoltaické simulaci proti metodě
Perturb & Observe (P&O). Prohledali jsme 20 000 přiřazení měřených signálů do složek oktonionu,
volby výstupní složky, znaménka a deformace metriky. Žádná konfigurace nepřekonala P&O. Ty nejlepší
se po symbolickém rozepsání redukují na rozhodovací pravidlo P&O, tedy „změna výkonu × změna směru“.
Jediná složka, kterou zveřejněný zdrojový kód skutečně počítá, `e0`, je pro libovolný vstup identicky
nulová; s nedeformovanou metrikou, kterou zveřejněný kód implikuje, se publikovaná konfigurace
umístila za 16 470 z 20 000 náhodně poskládaných variant. Dále ukazujeme, že asociátor tří po sobě
jdoucích stavů je rovný `[X, a, b]`, tedy trilineární funkci stavu a posledních dvou přírůstků, což
vysvětluje jak naměřenou třikrát pomalejší konvergenci, tak selhání při rychlých rampách ozáření:
užitečný signál se získává o jednu diferenci později a se o 74 % větším šumem měření. Veškerý kód je
k dispozici.

## 1. Co se tvrdí a co je zveřejněno

Preprint navrhuje dvouvrstvou architekturu: vrstvu „skryté kauzální paměti“ v oktonionové algebře
(O), která má nasčítávat tepelný drift, magnetickou saturaci a únavu materiálu, a observabilní
tenzorovou vrstvu tvořící střídu, přičemž most mezi nimi je oktonionový asociátor. Z tvrzených
důsledků jmenujme preemptivní ochranu hardwaru, vyřešení magnetické saturace a tepelné degradace
a bezdrátový fázový závěs roje převodníků přes společnou DC sběrnici.

Zveřejněné materiály obsahují:

| Co je potřeba k vyhodnocení návrhu | Je to tam? |
|---|---|
| pravidla násobení deformované algebry (všech 8 složek) | ne |
| pravidlo aktualizace deformace („metriky“) `g` | ne |
| fyzikální veličina, jednotka a normování každé složky `e0…e7` | jen popisky |
| definice argumentů X, Y, Z asociátoru | ne |
| převod asociátoru na střídu | v 1.0 jedna rozměrově nekonzistentní rovnice, v 1.1 nic |
| postup hledání MPP | ne |
| měření, grafy, porovnání s jakoukoli standardní metodou | ne |

Ve zdrojovém souboru `causal_associator.c` počítá funkce `OctonionMultiply` jen složku `e[0]`;
složky `e[1]` až `e[7]` jsou zastoupeny komentářem `// ...` a zůstávají neinicializované, takže
vrácený asociátor není číslo, ale nedefinované chování.

## 2. Metoda

**Algebra.** Implementovali jsme úplný oktonionový součin podle Fanovy roviny
(trojice (1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2), (7,1,3)) jako strukturní tenzor
a na 10⁵ náhodných trojicích jsme se strojovou přesností (maximální chyba ≤ 2,8·10⁻¹⁴) ověřili, že
implementace splňuje `e_i² = −1`, multiplikativitu normy `|xy| = |x||y|`, alternativitu
`[x,x,y] = [x,y,y] = 0`, flexibilitu `[x,y,x] = 0`, mocninnou asociativitu a že skalární složka
asociátoru je identicky nulová. Podle zveřejněného kódu se deformace aplikuje na každý součin po
složkách, `(ab)_k → g_k·(ab)_k`.

**Soustava.** Jednodiodový model 250 W 60článkového panelu (Isc 8,9 A, Voc 37,6 V, Rs 0,35 Ω,
Rsh 250 Ω), boost stupeň do 48 V baterie (`Vpv = Vbus(1−D)`), kvazistaticky v každém kroku regulace,
perioda MPPT 50 Hz, krok střídy 0,004 (asi 0,19 V), šum měření σ_V = 50 mV, σ_I = 20 mA. Skutečné
MPP se v každém vzorku hledá metodou zlatého řezu, takže účinnosti se vztahují ke skutečnému maximu,
ne k pevnému pracovnímu bodu. Použili jsme šest profilů ozáření a teploty (konstantní, skoky, rampy
ve stylu EN 50530, dva náhodné profily s mraky a teplotní drift okolí z 0 °C na 50 °C), dále model
částečného zastínění se třemi sekcemi s bypass diodami a přísný test se 100 ms hranami mraku
a rampami ±1000 W/m² během 2 s.

**Regulátory.** Oktonionový regulátor bere X, Y a Z jako tři po sobě jdoucí stavové oktoniony
(„kauzální paměť“) a střídu řídí z jedné složky `A_k` asociátoru, ve třech variantách: jen znaménko
(pevný krok, nejférovější srovnání s P&O), velikost kroku z `|A_k|` (1× až 10× základní krok, zesílení
vyladěné) a krok bez omezení. Referenční regulátory: P&O s pevným krokem, P&O s proměnným krokem
∝ |ΔP/ΔV| (zesílení laděné stejným postupem), P&O s periodickým projetím celé křivky, pevná střída
a náhodná procházka se stejnou velikostí kroku.

**Hledání.** 20 000 náhodných kandidátů nad: přiřazením 8 různých signálů z 12 (konstanta, V, I, P,
ΔV, ΔI, ΔP, T, Vbus, D, I/V, ΔD) do složek `e0…e7`; výstupní složkou; znaménkem; a buď `g = 1`, nebo
náhodnou teplotní deformací. Kandidáti byli hodnoceni na čtyřech trénovacích profilech a nejlepší
dvacítka byla přetestována na profilech, které neviděli.

## 3. Výsledky

### 3.1 Účinnost sledování MPP

Přísný test: průměr přes 6 profilů × 2 startovní body daleko od MPP (12 V a 33,6 V).

| Regulátor | s šumem | bez šumu |
|---|---|---|
| P&O, pevný krok | **99,35 %** | **99,45 %** |
| náhodná procházka, stejný krok | 50,9 % | – |
| konstantní napětí (střída pevně na 0,76·Voc) | 97,6 % | – |
| podle zveřejněného kódu, výstup `e0` | 49–94 % | 20–78 % |
| nejvěrnější rekonstrukce, výstup `e5` | 96,6 % | 97,2 % |
| jen jediný její člen, bez oktonionů | 97,1 % | 96,2 % |
| nejlepší z 20 000 prohledaných konfigurací | 98,4 % | 90,7 % |

Z 20 000 konfigurací překonalo 404 pevnou střídu a **žádná** nepřekonala P&O (medián 48,8 %,
99. percentil 98,9 %, maximum 99,78 % proti 99,81 % u P&O, hodnoceno při startu blízko MPP).
Třináct z nejlepší dvacítky ztrácí víc než 10 procentních bodů (nejvíc 61), když se odstraní šum
měření; několik z nich pak ujede na napětí naprázdno, kde je výstupní výkon nulový. Jejich sledování je tedy
částečně poháněné šumem.

Stojí za pozornost, jak dobrá je reference s konstantním napětím: regulátor, který nic neměří a jen
drží střídu na 0,76·Voc, dosáhne 97,6 %, tedy víc než nejvěrnější rekonstrukce návrhu (96,6 %) a víc
než všech 20 000 prohledaných konfigurací kromě 404 z nich. U jednoho panelu, jehož napětí MPP se přes
celý teplotní rozsah posune asi o 4 V, je taková reference obtížně překonatelná, a proto jako hlavní
referenci uvádíme P&O, ne tu naivní metodu.

**Citlivost na náhodné profily.** V lehkém pracovním bodě je odstup tak malý, že závisí na konkrétní
realizaci mraků. S jiným seedem generátoru u dvou mračných profilů překoná jedna konfigurace z 20 000
P&O o 0,02 procentního bodu při startu blízko MPP (99,80 % proti 99,79 %), zatímco v přísném testu ze
12 V a 33,6 V zůstává ta samá konfigurace 0,7 bodu za P&O s šumem a 8,5 bodu bez šumu. Stabilní je
tedy srovnání v přísném testu, které uvádíme; seedy generátoru u profilů jsou v kódu pevně zadané,
takže na nich žádné publikované číslo nezávisí.

Publikovaná konfigurace si zaslouží samostatný odstavec. S výstupem ze složky `e0` a nedeformovanou
metrikou, kterou zveřejněný kód implikuje (pravidlo aktualizace deformace není uvedeno nikde), dosáhne
za podmínek hledání 20,3 % účinnosti sledování, což ji řadí za 16 470 z 20 000 náhodně poskládaných
konfigurací (11,6. percentil). Naházet těch osm fyzikálních veličin do složek kostkou tedy vede
k lepšímu výsledku ve čtyřech případech z pěti. Doplnění libovolné teplotní deformace zvedne tu samou
složku na 97,0 %, pokud regulátor startuje blízko MPP, ale ze startu 12 V dá ten samý regulátor na
těch samých profilech jen 68,1 % a při konstantním ozáření 48,8 %, protože se od startovního bodu
prakticky nepohne. To zlepšení je tedy vlastností startovního bodu, ne metody.

### 3.2 Odkud berou funkční konfigurace svou informaci

U věrné rekonstrukce dá symbolické rozepsání složky `e5` osm bilineárních členů. Použije-li se jako
regulátor samotný člen `Vbus·(ΔG × Δ(ΔP))` (G = I/V), dosáhne 97,1 %, tedy o málo víc než celý
oktonionový výraz; ostatní členy přidávají jen šum. Nahrazení ΔP konstantou srazí regulátor na
33,9 %, nahrazení I/V na 25,0 %. U nejlepší konfigurace z hledání má rozpis podobu
`2·(Vbus·(ΔI × ΔP) − (ΔP × Δ(ΔD)))`, kde druhý člen je pravidlo P&O: změna výkonu násobená změnou
směru pohybu.

### 3.3 Proč konverguje pomaleji

Protože je asociátor trilineární a alternující, dávají tři po sobě jdoucí stavy
`[X, Y, Z] = [X, a, b]` s `a = Y − X` a `b = Z − Y` (ověřeno symbolicky), tedy funkci stavu
a posledních dvou přírůstků. Rozhodnutí tak stojí na druhých diferencích, jejichž šum je o 74 %
větší než u první diference (1,66 W proti 0,95 W v naší konfiguraci), a informativní část
determinantu se přitom z většiny vyruší, pokud se regulátor pohybuje konstantním krokem. Naměřený
podíl kroků učiněných správným směrem ve vzdálenosti větší než 2 V od MPP:

| Regulátor | správný směr | čistý postup za takt |
|---|---|---|
| P&O | **97,1 %** | 0,94 kroku |
| věrná rekonstrukce (e5) | 64,6 % | 0,29 kroku |
| nejlepší konfigurace z hledání | 68,7 % | 0,37 kroku |

Alternativita přidává strukturní nevýhodu: `[X, a, a] = 0`, takže výstup zaniká právě ve chvíli, kdy
regulátor vytrvale stoupá jedním směrem, a objeví se znovu teprve tehdy, když se trajektorie zakřiví.

### 3.4 Velikost kroku z „paměti“ a krok bez omezení

Odvození velikosti kroku z `|A_k|` (zesílení laděné přes deset řádů) dává dobu náběhu z 12 V 0,43 s
pro nejlepší konfiguraci proti 0,38 s u klasického P&O s proměnným krokem; ostatní konfigurace
potřebují 3,8 až 4,7 s. Úplné odstranění omezení kroku zhorší oktonionové regulátory z 98–99 % na
55–63 %, protože velikost asociátoru nemá jednotky ani kalibrovaný vztah ke vzdálenosti od MPP.
Ladění přidělilo P&O bez omezení nejmenší možný krok, tedy velké skoky zamítl sám ladicí postup.

### 3.5 Rychlé přechodové jevy a částečné zastínění

Na 100 ms hranách mraku nemá většina regulátorů co dohánět: napětí MPP se přes celý rozsah ozáření
posune asi o 1 V, zatímco 30 K teploty ho posune o 4 V. Výjimkou je věrná rekonstrukce, která na
rampě ±1000 W/m² za 2 s spadne na 72,6 % s průměrnou odchylkou napětí 8,55 V, protože změna ozáření
vnese do druhé diference velký falešný signál.

Při částečném zastínění (tři sekce, bypass diody; například 146 W při 17,6 V proti 80 W při 30,5 V
pro jednu zastíněnou sekci) nenašla globální vrchol žádná konfigurace, což se dalo očekávat: funkce
posledních dvou přírůstků nenese žádnou informaci o vrcholu vzdáleném 10 V. Pro porovnání, P&O
s periodickým projetím celé křivky našel globální vrchol ve 100 % úseků a ve stínu dosáhl 95,7 %
proti 87 % u obyčejného P&O.

### 3.6 Cena výpočtu

Jeden asociátor vyžaduje čtyři oktonionové součiny, tedy 256 operací násobení se sčítáním a osm
odečtení na krok regulace, proti jednomu porovnání u P&O.

## 4. Omezení

Vyhodnocujeme naši rekonstrukci, ne autorovu metodu, protože ta není definovaná; jiné přiřazení
veličin do složek se může chovat jinak, a právě v tom je problém, protože rozptyl mezi přiřazeními
je 20 % až 99 %. Model soustavy je kvazistatický a pokrývá jediný panel bez dynamiky převodníku;
dělení proudu mezi fáze, tvrzení o synchronizaci roje ani tvrzení o ochraně hardwaru netestujeme,
protože k nim neexistuje testovatelný popis. Naše výsledky neříkají nic o tom, jestli by oktonionová
vrstva mohla být užitečná v nějaké jiné formulaci.

## 5. Závěr

Pro samotné hledání MPP nepřináší oktonionový asociátor nic. Kde konfigurace funguje, funguje proto,
že rozpis náhodou obsahuje změnu výkonu, tedy veličinu, kterou P&O používá přímo, o jednu diferenci
dřív a s polovičním šumem. Složka, kterou zveřejněný kód počítá, je identicky nulová. Žádná
konfigurace z 20 000 nepřekonala metodu, jejíž princip je z konce 60. let.

Aby se návrh stal testovatelným, jsou potřeba tři věci: (i) tabulka násobení té algebry a pravidlo
aktualizace její deformace, s parametry, které jde identifikovat z měření; (ii) fyzikální veličina,
jednotka a normování každé složky, definice X, Y a Z a převod na střídu; (iii) jedno měřitelné
tvrzení porovnané se standardní metodou. Užitečný vlastní test pro jakýkoli navržený formalismus je,
jestli reprodukuje obyčejný případ: napiš v něm obyčejný převodník a ukaž, že z něj vyplynou
standardní rovnice.

Na závěr poznámka k metodě. Tvrzení vytvořená s pomocí velkého jazykového modelu působí sebevědomě
bez ohledu na to, jestli platí; model rozvine nepodloženou premisu stejně plynule jako správnou.
Nezávislé ověření, o jaké jsme se tady pokusili, je proto součástí práce, a ne dobrovolným doplňkem.
A je levné: celé tohle vyhodnocení je několik stovek řádků Pythonu. Formalismus navržený s takovým
modelem může být osmirozměrný, a přesto rozměrově nekonzistentní.

## Dostupnost dat a kódu

Veškerý kód, modely a experimentální skripty jsou k dispozici na https://github.com/karagos01/octonion-mppt-eval a reprodukují
každé číslo z tohoto komentáře (`run_all.sh`; hledání přes 20 000 konfigurací trvá asi 14 minut na
16 jádrech).

## Literatura

1. M. Mazgal, *Intensional Maxwellian Formalism in Power Electronics: Causal Octonion Control of
   Distributed Multi-Phase Converters*, Zenodo, 2026. DOI 10.5281/zenodo.22877912 (v 1.0),
   DOI 10.5281/zenodo.22914238 (v 1.1).
2. J. C. Baez, *The Octonions*, Bulletin of the American Mathematical Society 39(2), 2002, 145–205.
3. T. Esram, P. L. Chapman, *Comparison of Photovoltaic Array Maximum Power Point Tracking
   Techniques*, IEEE Transactions on Energy Conversion 22(2), 2007, 439–449.
4. N. Femia, G. Petrone, G. Spagnuolo, M. Vitelli, *Optimization of Perturb and Observe Maximum
   Power Point Tracking Method*, IEEE Transactions on Power Electronics 20(4), 2005, 963–973.
5. H. Patel, V. Agarwal, *Maximum Power Point Tracking Scheme for PV Systems Operating Under
   Partially Shaded Conditions*, IEEE Transactions on Industrial Electronics 55(4), 2008, 1689–1698.
6. D. P. Hohm, M. E. Ropp, *Comparative Study of Maximum Power Point Tracking Algorithms*,
   Progress in Photovoltaics: Research and Applications, 2003.
7. EN 50530, *Celková účinnost fotovoltaických střídačů připojených k síti* (dynamické testovací
   profily účinnosti MPPT).
