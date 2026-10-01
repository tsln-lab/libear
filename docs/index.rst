.. libear documentation master file

libear --- Recommendation ITU-R BS.2127 core library
====================================================

*libear* is a C++14 library to render ADM content according to :cite:`bs2127`.
It is not a complete application, but provides components to calculate gains
and apply them to audio, for embedding into applications which need to render
ADM content.

Rendering of ADM content is performed by two distinct processes:

- Calculating gains to apply to the input audio samples, as described in
  :doc:`calculating_gains`.
- Applying those gains to the input audio samples to produce output audio
  samples, as described in :doc:`dsp`.

To get started, check out the :doc:`installation` instructions.

Support
-------

DirectSpeakers, Objects and HOA typeDefinitions are supported, including
Cartesian (allocentric) positions, ``objectDivergence``, ``channelLock``,
``zoneExclusion``, ``screenRef`` and ``screenEdgeLock``. The following
parameters/features are currently not implemented:

HOA:

-  ``screenRef``
-  ``nfcRefDist``

Cartesian rendering requires the allocentric loudspeaker positions defined
for the ITU-R BS.2051 layouts; for other layouts only polar rendering is
available.

The Objects and DirectSpeakers gain calculators are checked against gains
produced by the reference implementation (the `EBU ADM Renderer
<https://github.com/ebu/ebu_adm_renderer>`_) over a large set of layouts and
block formats; see ``tools/reference/`` and
``tests/objects_reference_tests.cpp`` / ``tests/direct_speakers_reference_tests.cpp``.

*libear* does not include functionality to read BW64 files or parse ADM XML
data; for that functionality we recommend using libbw64_ and libadm_.

.. _libbw64: https://github.com/IRT-Open-Source/libbw64/
.. _libadm: https://github.com/IRT-Open-Source/libadm/

Credits
-------

*libear* is a joint development between IRT_ and `BBC R&D`_.

.. _IRT: https://lab.irt.de/
.. _`BBC R&D`: https://www.bbc.co.uk/rd

License
-------

.. code-block:: none

   Copyright 2019 The libear Authors

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

.. toctree::
  :maxdepth: 1
  :hidden:

  installation
  calculating_gains
  dsp
  loudspeaker_layouts
  conversion
  error_handling
  api/library_root
  changelog
  bibliography
