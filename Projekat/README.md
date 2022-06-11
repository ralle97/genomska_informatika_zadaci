# Burrows-Wheeler Transform and FM Index
Ovaj projekat predstavlja implementaciju osnosvnog Burrows-Wheeler algoritma i FMIndeksa, kao i daljih nadogradnji. Kod je pokriven sa unit i performance testovima. Takodje je napravljena i prezentacija koja analizira sve prethodno uradjeno kao i rezultate u vidu vremena izvrsavanja i memorijskog zauzeca izmedju osnovnog algoritma i nadogradjene varijante. 

# Pregled
* U data direktorijumu se nalaze 3 test seta nad kojima su sprovedeni benchmark testovi
* U memTest direktorijumu se nalaze testovi zauzeca memorije
* U sais direktorijumu se nalaze fajlovi za dalju nadogradnju algoritma pomocu SAIS SuffixArray
* U prezentacija direktorijumu se nalazi prezentacija implementacije algoritma, nadogradnji i rezultata
* Source code Burrows-Wheeler algoritma se nalazi u:
  * BurrowsWheelerTransform.ipynb
  * BurrowsWheelerTransformImproved.ipynb
  * BurrowsWheelerTransformSearchOverGenome.ipynb
* Unit testovi osnovnog Burrows-Wheeler algoritma se nalaze u BurrowsWheelerTransformTest.ipynb
* Benchmarking performance testovi se nalaze i pokrecu iz BurrowsWheelerTransformPerformaceTest.ipynb
* Osnovna implementacija FMindexa se nalazi u FMIndex.ipynb
* Nadogradnja FMIndexa kao i benchmarking i performance testovi se nalaze u FMIndexImproved.ipynb
* Dok se unit testovi nalaze u FMIndexTest.ipynb
* Kod za dalju nadogradnju pomocu SAIS algoritma se nalaze u SAIS.ipynb

# Instalacija
U jupyter notebook-u je potrebno pokrenuti !python -m pip install pydivsufsort, kako bi se instalirala biblioteka koja se koristi u kodu

# Video prezentacija
https://www.youtube.com/watch?v=hXuvTvzvfUE
