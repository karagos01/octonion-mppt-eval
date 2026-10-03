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
užitečný signál se získává o jednu diferenci později a se o 74 % větším šumem měření.

Dva výsledky na té rekonstrukci nezávisí vůbec. První: označení složek nenese žádnou informaci — z
5040 způsobů, jak přiřadit sedm fyzikálních veličin imaginárním jednotkám, nechá 168 násobicí tabulku
bitově identickou, bázová automorfní grupa má řád 1344 a těch 5040 přiřazení se zhroutí na 30
různých algeber, přičemž tatáž čtyři čísla vyjdou na dvou různých publikovaných orientacích. Čteno
s uvedenými jednotkami doslova staví 42 z 56 součinů do rovnosti neslučitelné dimenze a soustava
podmínek přes celou tabulku má hodnost 8 z 8, takže jediný dimenzionálně konzistentní oktonion nad ℝ
je ten, jehož složky jsou všechny bezrozměrné; `e0` je navíc multiplikativní jednotka, ne slot.
Druhý: *Maxwellian* z titulu je necitovaný a atribuce sleduje algebru, ne historii — ze všech deseti
autorových depositů na Zenodu je Maxwell jmenovaný jen ve čtyřech, které ho párují s oktoniony,
nulakrát ve třech kvaternionových, a „Treatise", „1873", „1865" a „Heaviside" se nevyskytují nikde.
Veškerý kód je k dispozici.

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

## 4. Označení složek nenese žádnou informaci

Sekce 2.1 ve v1.1 je jediné místo, kde do oktonionové vrstvy vstupuje fyzika, takže celý návrh stojí
na jedné větě:

> „Let O be a state octonion mapping dimensions e_0 … e_7 to physical potentials (e.g., scalar
> voltage, vector current, thermodynamic entropy, and electromotive gradients)."

`labels.py` měří na té větě tři věci. Žádná z nich nezávisí na naší rekonstrukci regulátoru, takže
tato sekce stojí nezávisle na sekci 3.

**Těch sedm imaginárních jednotek je záměnných.** Fyzikální přiřazení je bijekce ze sedmi veličin na
e_1 … e_7, takže jich je 7! = 5040. Dvě přiřazení popisují tutéž algebru, kdykoli je přeznačení mezi
nimi automorfismus násobicí tabulky. Spočítáno přímo, na orientaci používané v celém tomto
repozitáři i na té, kterou tentýž autor později publikoval jako SOTP Eq. (3):

| | naměřeno |
|---|---|
| přeznačení, po kterých je tabulka bitově identická | **21 z 5040** |
| …když se smí obrátit znaménko každé jednotky | **168 z 5040** |
| řád celé bázové automorfní grupy | **1344** = 8 × 168 |
| různých algeber dosažitelných přeznačením | **30** |

Obě tabulky dávají tatáž čtyři čísla, takže je to vlastnost oktonionů, ne jedné konkrétní orientace.
Každé z těch 5040 přiřazení je algebraicky identické se 167 ostatními. Hlubší důvod je, že Aut(𝕆) nad
reálnými čísly je čtrnáctirozměrná výjimečná grupa G₂, která působí tranzitivně na bázových trojicích:
žádný invariant algebry nerozliší e_1 od e_4. Pojmenovat jednu složku „thermodynamic entropy" a jinou
„electromotive gradient" tedy nepřidává žádnou podmínku, kterou by jakýkoli výpočet mohl vidět. Je to
v souladu s výsledkem našeho prohledávání v sekci 3.2 — tam, kde konfigurace fungovala, fungovala
kvůli tomu, *který signál* se přivedl, nikdy kvůli tomu, *do kterého slotu*.

**Žádné přiřazení jednotek není možné.** Čteno doslova, s e_0 ve voltech, e_1 … e_3 v ampérech, e_4
v J/K a e_5 … e_7 ve V/m, **42 z 56** součinů dvou různých imaginárních jednotek staví do rovnosti
neslučitelné dimenze: e_1 e_2 = e_4 vyžaduje [A][A] = [J/K], e_1 e_3 = e_7 vyžaduje [A][A] = [V/m] a
tak dál. To jsou všechny.

Silnější tvrzení je, že to nezachrání žádné přeznačení ani přeškálování. Napiš jednu neznámou
dimenzionální exponentu d_i na složku a seber podmínku d_i + d_j = d_k ze všech 64 položek tabulky.
Výsledná soustava má **hodnost 8 v 8 neznámých**, takže její prostor řešení má dimenzi **0**: jediný
dimenzionálně konzistentní oktonion nad ℝ je ten, ve kterém je všech osm složek bezrozměrných. Proč,
se vidí na dvou položkách — e_0 e_i = e_i vynutí d_0 = 0 a e_i e_i = −e_0 pak vynutí d_i = 0. Oktonion
je algebra nad ℝ; jeho osm složek se k sobě sčítá, takže musí nést jednu dimenzi. Sekce 2.1 nepopisuje
bohatší objekt než vektor osmi reálných čísel. Popisuje výraz, který nelze vyhodnotit.

**e_0 není slot.** e_0 e_0 = e_0: je to multiplikativní jednotka. Přiřadit jí „scalar voltage" tvrdí,
že V² = V. Mapování nabízí osm jmen pro sedm záměnných míst a jedno místo, které místem není.

Toto je ta dimenzionální nekonzistence, na kterou sekce 1 odkazuje jako na „jednu dimenzionálně
nekonzistentní rovnici" ve v1.0, nalezená u zdroje, ne v odvozeném vzorci.

## 5. Jaká algebra se připisuje Maxwellovi

Titul obou verzí je *Intensional Maxwellian Formalism* a abstrakt v1.0 popisuje práci jako „restoring
James Clerk Maxwell's original hypercomplex architecture". To tvrzení je nosné: dodává historické
oprávnění, proč vůbec volit nonasociativní algebru. A je necitované — preprint neodkazuje na žádnou
Maxwellovu práci a nemá bibliografii vůbec.

`census.sh` stáhne všech deset autorových depositů ze Zenoda, převede každý na text a spočítá klíčová
slova. Je to jediný skript tady, který potřebuje síť; počty níže jsou z běhu 3. října 2026.

| deposit | konec DOI | Maxwell | octonion | quaternion | Treatise | 1873 | 1865 | Heaviside |
|---|---|---|---|---|---|---|---|---|
| akustický projektor | 22876757 | 2 | 19 | 0 | 0 | 0 | 0 | 0 |
| X-Ternary | 22877345 | 1 | 2 | 0 | 0 | 0 | 0 | 0 |
| MPPT v1.0 | 22877912 | 4 | 13 | 0 | 0 | 0 | 0 | 0 |
| MPPT v1.1 | 22914238 | 5 | 12 | 0 | 0 | 0 | 0 | 0 |
| COBAR | 22923154 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| topologický pohon | 22934058 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CQFT | 22959886 | **0** | 0 | 10 | 0 | 0 | 0 | 0 |
| PCTP | 22962188 | **0** | 0 | 15 | 0 | 0 | 0 | 0 |
| SOTP | 22962557 | **0** | 20 | 8 | 0 | 0 | 0 | 0 |
| openQL / openOL | 23112560 | 2 | 3 | 7 | 0 | 0 | 0 | 0 |

Z tabulky vypadnou tři vzorce.

**Maxwell je jmenovaný jen tam, kde jsou oktoniony.** Všechny čtyři deposity, které Maxwella párují
s inženýrskou aplikací, mu připisují dvouvrstvou architekturu, jejíž první vrstva je oktonionová:
„Following Maxwell's formalism, the system's underlying physical models are structured into two
distinct layers. The first layer is the Hidden Causal Memory, utilizing Octonion Algebra" (akustický
projektor); „a two-layer formalism inspired by Maxwell's original equations … Layer 1 (Hidden Causal
Memory): The system's causal history and interaction sequence are stored using non-associative
octonion algebra" (X-Ternary); „restoring James Clerk Maxwell's original hypercomplex architecture"
(MPPT v1.0).

**Kde je algebra kvaternionová, Maxwell zmizí.** CQFT a PCTP jsou postavené na kvaternionech a
jmenují ho nulakrát. SOTP, oktonionový člen téže trilogie, ho jmenuje taky nulakrát. Místo
historického oprávnění nasazuje CQFT jiné — Adlerovu kvaternionovou kvantovou mechaniku — a tvrdí, že
„this historical failure was not due to the nature of quaternions, but rather the flawed axiom of
scalar fungibility". Autorita, která se připne přesně tam, kde jsou oktoniony, a odpadne, když se
algebra změní, je vybíraná k závěru, ne konzultovaná.

**Necituje se nic.** „Treatise", „1873", „1865" a „Heaviside" se ve všech deseti depositech vyskytují
nulakrát. Ať Maxwell původně napsal cokoli, žádný deposit neříká kde.

To, co historický záznam obsahuje, je užší než to tvrzení a míří na tu druhou algebru. *Dynamical
Theory of the Electromagnetic Field* z roku 1865 je dvacet skalárních složkových rovnic ve dvaceti
proměnných, bez hyperkomplexní algebry jakéhokoli druhu. Kvaterniony se v *Treatise* z roku 1873
objevují v §§618–619, v posledních dvou článcích kapitoly IX, pod záhlavím „Quaternion Expressions
for the Electromagnetic Equations" — na necelých dvou stranách dvousvazkového díla o přibližně tisíci.
§618 je jmenný seznam vektorů a skalárů, které už byly implicitně v užívání; §619 říká, že *pokud* se
„vector" a „scalar" čtou jako vektorová a skalární část kvaternionu a *pokud* se ∇ bere jako
kvaternionové, pak se už odvozené rovnice dají zapsat jako (A)–(L). Vyskytují se tam jen operátory
`S.` a `V.`. Maxwell nikde v celé knize nevynásobí dva kvaterniony: plný součin PQ se neobjevuje. Je
to Hamiltonova operátorová notace, ne kvaternionová algebra, a reálná část kvaternionu nenese žádnou
fyziku.

Maxwellova vlastní pozice je doložená. „Endeavoured to avoid any process demanding from the reader a
knowledge of the Calculus of Quaternions", doporučoval „the introduction of the ideas, as
distinguished from the operations and methods of Quaternions" a Taitovi napsal, že chce „leaven my
book with Hamiltonian ideas without casting the operations into a Hamiltonian form". Měl i fyzikální
námitku: vektorový kvaternion dá po umocnění minus kvadrát své délky, takže kinetická energie zapsaná
kvaternionově vyjde negativní. V kapitole X se gotické symboly objevují zbavené kvaternionového čtení
a potom už vůbec.

Takže silná forma toho tvrzení není dostupná ani v jedné algebře. Oktoniony — Graves 1843, Cayley
1845 — se u Maxwella neobjevují nikde, ani jako notační poznámka. A nic potlačené nebylo: ty dvě
stránky jsou pravděpodobně nejvlivnější v celé knize, protože z nich Heaviside a Gibbs postavili
moderní vektorovou analýzu, a sám Maxwell zavedl „gradient", „convergence" a „curl" ve své eseji
o klasifikaci fyzikálních veličin z roku 1871. Co z kvaternionové formy vypadlo, byla reálná část,
kterou Maxwell nikdy nenaplnil fyzikou.

Pro úplnost: openQL / openOL (3. října 2026) je první deposit, který Maxwella připíná ke
kvaternionům, a ne k oktonionům — „4D quaternion mechanics (for computing phase shifts, rotations,
and Maxwell's equations)". To je ta algebra, u které má historické tvrzení aspoň slabý základ. Pořád
je necitované a pořád bez mechanismu, ale atribuce teď míří na správnou algebru. Ten deposit je
vyhodnocený samostatně v `openql-eval`.

## 6. Omezení

Vyhodnocujeme naši rekonstrukci, ne autorovu metodu, protože ta není definovaná; jiné přiřazení
veličin do složek se může chovat jinak, a právě v tom je problém, protože rozptyl mezi přiřazeními
je 20 % až 99 %. Model soustavy je kvazistatický a pokrývá jediný panel bez dynamiky převodníku;
dělení proudu mezi fáze, tvrzení o synchronizaci roje ani tvrzení o ochraně hardwaru netestujeme,
protože k nim neexistuje testovatelný popis. Naše výsledky neříkají nic o tom, jestli by oktonionová
vrstva mohla být užitečná v nějaké jiné formulaci.

## 7. Závěr

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
16 jádrech).  Sekce 4 a 5 jsou ty levné: `labels.py` reprodukuje všechny počty
ze sekce 4 asi za sekundu a nepotřebuje nic než numpy a `census.sh` reprodukuje tabulku ze sekce 5
stažením těch deseti depositů ze Zenoda — je to jediný skript tady, který používá síť.

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
8. J. C. Maxwell, *A Treatise on Electricity and Magnetism*, Clarendon Press, 1873, §§618–619
   („Quaternion Expressions for the Electromagnetic Equations").
9. J. C. Maxwell, *A Dynamical Theory of the Electromagnetic Field*, Philosophical Transactions of
   the Royal Society 155, 1865, 459–512.
10. J. C. Maxwell, *On the Mathematical Classification of Physical Quantities*, Proceedings of the
   London Mathematical Society 3, 1871, 224–233 (kde zavádí „gradient", „convergence" a „curl").
11. N. Wheeler, *Theories of Maxwellian Design*, Reed College, 1998 (k §§618–619, korespondenci
   Maxwell–Tait a námitce o negativní kinetické energii).
12. J. M. Chappell, A. Iqbal, J. G. Hartnett, D. Abbott, *The Vector Algebra War: A Historical
   Perspective*, IEEE Access 4, 2016, 1997–2004. arXiv:1509.00501.
13. A. Hurwitz, *Über die Composition der quadratischen Formen von beliebig vielen Variabeln*,
   Nachrichten der Gesellschaft der Wissenschaften zu Göttingen, 1898, 309–316.
