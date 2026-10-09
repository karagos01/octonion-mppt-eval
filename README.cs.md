# Oktoniony v MPPT — nezávislé numerické vyhodnocení

Reprodukovatelné vyhodnocení tvrzení, že oktonionový asociátor `[X,Y,Z] = (XY)Z − X(YZ)` umí řídit
fotovoltaický MPPT převodník, jak navrhuje M. Mazgal, *Intensional Maxwellian Formalism in Power
Electronics* (Zenodo, DOI 10.5281/zenodo.22877912 a 10.5281/zenodo.22914238) a repozitář
`Causal-Octonion-MPPT.`

Text je v [`PAPER.cs.md`](PAPER.cs.md), anglicky v [`PAPER.md`](PAPER.md). Krátce: kompletní
oktonionový součin je implementovaný a ověřený, chybějící definice jsou rekonstruované co nejvěrněji,
prohledáno je dalších 20 000 konfigurací a žádná nepřekonala Perturb & Observe. Složka, kterou
původní kód počítá (`e0`), je pro libovolný vstup identicky nulová.

> Komentáře v kódu a výpisy skriptů jsou anglicky. Anglická verze textu je v [`README.md`](README.md)
> a [`PAPER.md`](PAPER.md).

## Co je potřeba

```
python3 -m pip install -r requirements.txt   # numpy, sympy
```

GPU není potřeba. `search.py` používá 16 procesů (`multiprocessing.Pool(16)`), na menším počtu jader
tu hodnotu sniž.

## Reprodukce výsledků

```
./run_all.sh          # všechno, asi 20 minut na 16 jádrech
```

nebo jednotlivě, v tomto pořadí (pozdější skripty čtou `search.npz`):

| Skript | Co spočítá | Čas |
|---|---|---|
| `octonion.py` | ověření algebry (norma, alternativita, skalární složka asociátoru) | 3 s |
| `labels.py` | jestli označení složek ze sekce 2.1 nese informaci: bázové automorfismy, dimenzionální konzistence | 1 s |
| `baseline.py` | P&O, pevná střída, náhodná procházka a všech 16 věrných variant na 6 profilech | 4 min |
| `symbolic.py` | důkaz, že `[X,Y,Z] = [X,a,b]`, a rozpis každé složky symbolicky | 20 s |
| `convergence.py` | konvergence z 12 V a z 33,6 V, s šumem měření i bez něj | 1 min |
| `policy.py` | rozhodovací tabulka bez šumu: ví regulátor, kde je MPP? | 30 s |
| `search.py` | hledání přes 20 000 konfigurací, zapíše `search.npz` | 14 min |
| `percentile.py` | kde se umístila publikovaná konfigurace mezi 20 000 náhodnými | 30 s |
| `strict.py` | přísné přetestování nejlepší dvacítky, zastoupení signálů u úspěšných | 2 min |
| `ablation.py` | který signál a který jediný člen nese informaci | 3 min |
| `adaptive.py` | velikost kroku z velikosti asociátoru, ladění zesílení | 1 min |
| `adaptive2.py` | doba náběhu, ustálená účinnost, zastínění přes 8 průběhů šumu | 2 min |
| `free.py` | krok bez omezení; reference P&O s periodickým projetím křivky | 3 min |
| `clouds.py` | 100 ms hrany mraku a rampy ±1000 W/m² za 2 s | 2 min |
| `why_slow.py` | podíl kroků správným směrem; šum v první a druhé diferenci | 1 min |

## Moduly

- `octonion.py` — strukturní tenzor Fanovy roviny, součin, asociátor, deformace `g` po složkách,
  vlastní testy.
- `sim.py` — jednodiodový model panelu, boost stupeň, profily ozáření a teploty, regulátory
  (P&O, P&O s proměnným krokem, P&O s projetím křivky, pevná střída, náhodná procházka, oktonionový
  asociátor v režimech znaménko / proměnný krok / bez omezení), simulační smyčka vektorizovaná přes
  N regulátorů.
- `shading.py` — tři sekce s bypass diodami, křivka P(V) s více vrcholy, globální MPP.
- `census.sh` — **jediný skript, který potřebuje síť.** Stáhne všech deset autorových depositů ze
  Zenoda, převede je přes `pdftotext` a spočítá klíčová slova, na kterých stojí atribuce Maxwellovi
  (PAPER.cs.md §5).

## Klíčová čísla

| Regulátor | s šumem | bez šumu |
|---|---|---|
| P&O, pevný krok | **99,35 %** | **99,45 %** |
| náhodná procházka | 50,9 % | – |
| konstantní napětí (střída na 0,76·Voc) | 97,6 % | – |
| podle zveřejněného kódu, výstup `e0` | 49–94 % | 20–78 % |
| nejlepší věrná rekonstrukce (`e5`) | 96,6 % | 97,2 % |
| jen jediný její člen, bez oktonionů | 97,1 % | 96,2 % |
| nejlepší z 20 000 konfigurací | 98,4 % | 90,7 % |

Průměr přes 6 profilů ozáření a teploty × 2 startovní body daleko od MPP. Z 20 000 konfigurací
překonalo 404 pevnou střídu a žádná nepřekonala P&O.

Publikovaná konfigurace (výstup `e0`, nedeformovaná metrika, jak plyne ze zveřejněného kódu) dosáhne
20,3 % a umístila se za 16 470 z 20 000 náhodně poskládaných konfigurací. Naházet těch osm veličin do
složek kostkou ji tedy překoná ve čtyřech případech z pěti (`percentile.py`).

Kostka si vede stejně dobře proto, že ty sloty jsou záměnné. Ze 7! = 5040 způsobů, jak namapovat sedm
fyzikálních veličin na `e1…e7`, nechá **168** násobicí tabulku bitově identickou, takže každé
přiřazení je algebraicky nerozlišitelné od 167 dalších a těch 5040 se zhroutí na **30** různých
algeber; celá bázová automorfní grupa má řád **1344** (`labels.py`, tatáž čtyři čísla na orientaci
z tohoto repozitáře i na té, kterou autor později publikoval jako SOTP Eq. 3). Čteno s jednotkami ze
sekce 2.1 doslova staví **42 z 56** součinů různých imaginárních jednotek do rovnosti neslučitelné
dimenze a soustava podmínek `d_i + d_j = d_k` přes celou tabulku má hodnost **8 z 8** — takže jediný
dimenzionálně konzistentní oktonion nad ℝ je ten, jehož všech osm složek je bezrozměrných. `e0` navíc
není slot vůbec: je to multiplikativní jednotka, takže nazvat ji „scalar voltage" tvrdí, že V² = V.

*Maxwellian* v titulu je necitovaný a `census.sh` ukazuje, že atribuce sleduje algebru, ne historii:
Maxwell je jmenovaný ve čtyřech depositech, které ho párují s oktoniony (2, 1, 4 a 5krát), a
**nulakrát** v CQFT, PCTP a SOTP, kde je algebra kvaternionová. „Treatise", „1873", „1865" a
„Heaviside" se ve všech deseti depositech vyskytují **nulakrát**. Kvaterniony jsou v *Treatise* z roku
1873 v §§618–619 — na dvou stranách hamiltonovské operátorové notace, ve kterých Maxwell nikdy
nevynásobí dva kvaterniony — a oktoniony se u Maxwella neobjevují vůbec nikde.

## Doprovodná vyhodnocení

- [`xternary-eval`](https://github.com/karagos01/xternary-eval) — 2bitový inferenční engine pro LLM
- [`causal-trilogy-eval`](https://github.com/karagos01/causal-trilogy-eval) — trilogie CQFT / PCTP / SOTP ze září 2026
- [`openql-eval`](https://github.com/karagos01/openql-eval) — tenzorová maticová architektura openQL / openOL z října 2026
- [`cymatic-eval`](https://github.com/karagos01/cymatic-eval) — cymatická stimulace průduchů a pulzní osvětlení z října 2026
- [`reference-audit`](https://github.com/karagos01/reference-audit) — jestli všech 69 citací ve všech 16 depositech říká to, pro co je citovaných

## Licence

Kód (všechny `*.py` a `run_all.sh`): MIT, viz `LICENSE`.
Text `PAPER.md` a `PAPER.cs.md`: CC BY 4.0.
