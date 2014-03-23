WORLD - a high-quality speech analysis, modification and synthesis system

WORLD is free software for high-quality speech analysis, modification and synthesis.
It can estimate Fundamental frequency (F0), spectral envelope and excitation signal,
and also generate the sound like input voice with only estimated parameters.

2. Usage
Please see test.cpp.

(1) F0 estimation by Dio()
(1-1) F0 is refined by StoneMask() if you need more accurate result.
(2) Spectral envelope estimation by Star()
(3-1) Excitation signal extraction by Platinum()
(3-2) You can use AperiodicityRatio() if you want to use 
      aperiodicity instead of extitation signal.
(4) You can manipulation these parameters in this phase. 
(5-1) Voice synthesis by Synthesis(), provided that Platinum() is used.
(5-2) Voice synthesis by SynthesisFromAperiodicity(), provided that AperiodicityRatio() is used.

English document is written by a Japanese poor editor.
I willingly accept your indication on my English text.

3. License
WORLD is free software, and you can redistribute it and 
modify it under the terms of the BSD License.
Please see copying.txt for more information.
You can use this program for business, while I hope that 
you contact me after you developed software with WORLD.
This information is crucial to obtain a grant to develop WORLD.

4. Contacts
WORLD was written by Masanori Morise.
You can contact him via e-mail (mmorise@yamanashi.ac.jp)
or Twitter: @m_morise

