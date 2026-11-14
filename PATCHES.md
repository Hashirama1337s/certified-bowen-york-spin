# PATCHES: lab paths in the shipped copies

The lab's files name lab folders in a few places. To keep machine paths out of the public repository, every
shipped copy has had those prefixes replaced by placeholders from one fixed table, `tools/path_patch.json`,
applied longest-first. Nothing else in any file was changed. The substitution is exactly reversible:
`tools/check_pins.py` restores every patched file in memory and checks that the restored bytes have the original
sha256 (for sealed files: the lab ledger pin), and that the shipped file has the patched sha256 listed here.
`tools/reproduce.py` restores the originals inside its work folder before running anything, so every run sees the
sealed bytes and every pinned input passes its sha256 check. The original strings are stored in the table as hex
so that the table itself carries no path text.

## Substitution table (applied in this order)

| Entry | Placeholder | What the original names |
|---|---|---|
| E1 | `<LAB>\\` | lab project root, JSON-escaped Windows form (inside recorded messages) (11 characters) |
| E2 | `<SCRATCH>` | the word for the lead's temporary working folder (in disclosures and verdict notes) (10 characters) |
| E3 | `<LAB>/` | lab project root, forward-slash form (docstrings, usage lines, engine-location literal, test fixtures) (9 characters) |
| E4 | `<lab>/` | lab project root, lower-case form (scratch ledger rows inside self-sha controls) (9 characters) |

## Patched files

| File | Pinned | Original sha256 | Patched sha256 | Entries (count) |
|---|---|---|---|---|
| code/b3t.py | yes | `0043585b666fd0f5b0afd2ab4f16dc98335600834943c6a8bd7fdee13aed679a` | `74cd23ae2ac5f1925c8fa08518fb17c49e6ea86428cc2ab08640dbef32169211` | E3 x2, E4 x4 |
| code/fluxt.py | yes | `957929cc1d4f3111f1f3d4ec33ffb0df7bed5e94f8465a9aca03871a10fa0c27` | `4af8168eead8cbc58a39b1bfd18c2f7a4a0c8ea519cb40bfb1cd14223cac2cd4` | E3 x2 |
| code/gravity_engine/bowen_york_sqrt.py | yes | `6f3a40993bf6ca4b082ee7699913fd49dec3beecf15365dfa12e730ccc7bc9b3` | `fdad2ec01e79846dd71adda71c34af618c8d800dc368290ba97bfdbc3f95f109` | E3 x6 |
| code/gravity_engine/lichnerowicz_axisym.py | yes | `3659ae21df8534a74ede3af11a1d7fbae670b7478d67c0bb2bd5b34710d79f52` | `a55e17db1d8eb6663e811d5599d2eb2338eafe112bc188a4eb327f92ef448d26` | E3 x3 |
| code/gravity_engine/tests/_v2_node_child.py | no | `296f1717211f2db4f53facfb6adf811218316e1df92c9036ee009e89888e5155` | `ae9d03db8601f79c5282c00c25c874a511f6b9829ba48d21acb2b0c0e49aeb5d` | E3 x1 |
| code/gravity_engine/tests/test_bowen_york_sqrt.py | no | `5f274b673c6814663e55d92faa60646578e0b3df5fe5c9d4ffb6d83a3f63b186` | `1fa654ae3bc2f7144b6da68c5ae52687a371db412d6b405890a591b98fea5cf9` | E3 x7 |
| code/gravity_engine/tests/test_bowen_york_tangent.py | no | `fccd630bb109808e720828b92afc4a1ceb3d58a09a211efa0f3ec8f1aed0681f` | `432029f4d8658c81a14e5bbebc81eb15f4e3f48ec74bd449f133ab523c70035a` | E3 x4 |
| code/gravity_engine/tests/test_engine_v2.py | no | `233fe7e26fb69bc9dd8bf349e5defcab23eb5378507d45b5cc397724039a91d6` | `2de4b7031b4a9c351c35f772a0defde100f8f56c7cc400764a9f3839453dcddc` | E3 x4 |
| code/gravity_engine/tests/test_lichnerowicz_axisym.py | no | `35bf10072649606cc5670f5f0c990f9883d692b05f6e9134f502bf22ab25287c` | `cf24c59e29ae36a2b6afe82aa221bbb1982be3ed21c7e566415c9fd7f7a5576e` | E3 x6 |
| code/p1.py | yes | `23c6cbfdae63d819bbb75a5da75bece0d1c6e65c77ca36846448652626fa810f` | `af94ebf0536045b09719dd95915f25473629d7f0cc40eb1a56bc87bd7a87cf74` | E3 x2 |
| code/p2a.py | yes | `feee7d6f273fbfe97581dca0d224c0e387b6cdf31f9f10ae6f949543de2dbf04` | `cfe84510e300b1ad623f50e003ccce3af45456814c473ec191d0771f742e4b74` | E3 x2 |
| code/p2b.py | yes | `124ce51a261935c5cdc80f9956712d3700cefd0d1177712d09d944a05afa7f2d` | `48761b41c4fcb9d9fc51928f74a70ed972e783c4650acbe6f20e0a465f7e0d2b` | E3 x4 |
| code/p2c.py | yes | `9d2452bd865bffbe1742454e1e19dab2882d10ffdc7dbd83d5e0ade985fd311b` | `fb698140fbd04120944877e565af6d43401fa74296ef5b261fda2789c96ac3f3` | E3 x4 |
| code/p2d.py | yes | `f1fae7330e2472ac4f5e558ad40c88ea03ac15d052181c2952ef307ec1e75a1b` | `539758f63e2676c80db06e513038ba46d6aac881da40425554490e2b4b1fe8d5` | E3 x3, E4 x2 |
| code/p2e.py | yes | `a25e6cc4fc687bfb93e3a354a4010e4e6bd71de8de6aaeac57229ac7b4873555` | `434b7b5aabef4ac43c55f17f4f0f737b2b85bf4d9e37bb66169947e9ee2b01ad` | E3 x3, E4 x3 |
| code/rep/rep5/rep5.py | yes | `9b003febfedd1c5765f01f1ddff10a8b01c11f7550125b9adf9a6c787f02cb89` | `17a7d518e97a8d1c14657a7198a78d24513f09fd3728cbefbc44a0834502a40d` | E3 x2 |
| code/rep/rep7/rep7.py | yes | `c473c57be1f92f963af66b4d27e236d1c3f86d9541d324c1b36af41eaa019c32` | `c871100e65cb70fef73aef4d3794e96cf1841f69472ea6ef091949486602f8d2` | E3 x1 |
| results/b3t/controls_b3t.json | yes | `cab8efe8e828ce9c4c7b2b5f832ca7627d97cc21df1b76c1da0474abe6ff683a` | `66da26be04a98e89513bb32dcc3a3ef628e3bd378ab04ef60e3952c7341026d1` | E1 x1 |
| results/ledger-bhslice.json | no | `e7130980249d6ddcbade4411088b8878199c76768065fbe50be852409c59580b` | `11ec29ccc2b3d79b75119fb9a38e5eee841d906e7170dce3a06ded620d808e1f` | E2 x1 |
| results/p2a/TARGETS-P2.md | yes | `e42e8b1e0adf37e9e1ae8dd9d3a2362a2b349584072877326ad24edd6233cbce` | `436d0e0653d568e62955d80e124081ebf4a6627007daf5cfd024d62541eb37ca` | E3 x1 |
| results/p2a/controls_p2a.json | yes | `481e7a8c2667318b4adfd9e87e6c93a722906c848dadaae6df312e76d66e3a49` | `5ae6cc62ea5b453ed49e1364993a94a85203a9804fd3aa2bd540084fda39fad5` | E3 x1 |
| results/p2b/controls_p2b.json | yes | `2508deae71c159b07c6b48d43ff72ecccdef7baf4796d1a045ea9685cf85fe35` | `651aca37561a7a4ffada1bbda100e00103151990d12ea73ff86648db3bc6993e` | E3 x1 |
| results/p2b/result_p2b.json | yes | `9f4fa988dfb66bb6978f7e41fe6b26d32c0c152b908061ffe91076bc154776a3` | `762146d0c532c58b2db83507c8a40c0736e9594641fdd589bd26aeec6ac321cd` | E3 x1 |
| results/p2c/controls_p2c.json | yes | `a87d42f1d60b9c62fca504ffe86d54913d64f9e48fa8095ccaa72797bde29518` | `5a5153b44fe10516c19b733c5e09a42e832329b61ab6718955c27e5cf815cd1e` | E1 x1, E3 x1 |
| results/p2d/AMEND-P2d-3.md | yes | `3afd48931b58d95a28b73c140be06efc65235c370eaa1133f9714ab68e3b3850` | `861ac0ee18876f4409994da816bc824008a39f0697250affa3cc63a1eeb41dc3` | E2 x1 |
| results/p2d/TARGETS-P2d.md | yes | `5f76f324c361c3c92e557927f28c487439cf44756916c44f5257a619d88d4865` | `c5015594ee2a140b7d158687baeba4d952aa5b3e2dce0032587383dd5df91e8e` | E3 x1 |
| results/p2d/VERDICT-P2d.md | no | `e42a1fd817aa5cb5a9cdb3b829c1ddda1e2196ed93d8c2f8a27a83f1b45347f7` | `56d50ed5e762b72c3acfc51e6b9c2cc3a22ac616c0ed55e806eba7cd0705d652` | E2 x2 |
| results/p2d/controls_p2d.json | yes | `9a638fda7c84c430732146b5cbce50393366af5208457753d75d4cd019e2d045` | `6daad7cee6b94178f30fe29870043136eea1d407970c31569001f47d751497c3` | E1 x1 |
| results/p2e/VERDICT-P2e.md | no | `ea988780f25e6b553ae86600826e4cc7223fd7455bb1aa5c6aad123314f744b0` | `8d8e5edf68804564f75ee69d02048b3966860889451b98fe576c75a30e8498d0` | E2 x1 |
| results/p2e/controls_p2e.json | yes | `e96bc6295b0e1a942ef60d1d5abaf2c08fa2252a22bbee5b6ef6ff77897ff474` | `e8a05bedc6c8d7c90e57f9c507b44ec251b4bacb5fc962a963575aae536508e8` | E1 x1 |
| results/rep/rep5/RESULT.md | no | `249f9dfce34e728ad73b66b0c4ae9a58637776487e6488618d42ce628adeb966` | `c019963b36cf01c51eb90d3928b1e79191037a516a9831d07f1cbf7bb1619804` | E2 x2, E3 x9 |
| results/rep/rep5/result_rep5.json | no | `593db63eb381b5160159833fb805e066e1db5d9c230710a11a654cf67f4b8358` | `67fb429429c7038ae738e41a0b2fd9a01f274c39f9664948bf5ad7c9b18e06d5` | E3 x1 |

Pinned = the original sha256 is pinned by the lab ledger or by a later instrument (`tools/pins.json`). For the
other files the original sha256 is the lab file as it stood when this repository was assembled; it is listed so
that the restoration can be checked the same way. No file was redacted: every change above is reversible.

Separately (not a path patch): `results/b2v/b2v_result.json` is the JSON object of the lab file without the shell
residue a redirect appended (`results/b2v/NOTE.md`).
