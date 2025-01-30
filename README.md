# FragmentTools

## 概要
ABINIT-MP による FMO 計算で、手動フラグメント分割をサポートするツール群


## ツール群
* `fred4.py`
	: .ajf → .fred 間、.fred → .fred 間、.fred → .ajf 間でファイルを相互変換するプログラム。
	* `autofrag`
		: .pdb ファイルを基にフラグメント分割をして、.fred ファイルに出力する。
	* `edit`
		: .ajf → .fred の変換をする。
	* `rewrite`
		: .fred → .fred の変換をする (フラグメント編集後に番号を振り直す)。
	* `output`
		: .fred → .ajf の変換をする。
	* `autofrag`
		: `autofrag.py` の機能
	* `editfrag`
		: `editfrag.py` の機能
* `editfrag.py`
	: フラグメント情報を構造を与えてフラグメントを編集する。
	: 複数の分子に対して同じフラグメント分割をする。


## 使用方法
### fred4.py
#### autofrag mode
```sh
$ fred4.py autofrag [-h] -p INPUT.pdb -o OUTPUT.fred [-sp FRAGMENT_OPTION_PROTEIN] [-sn FRAGMENT_OPTION_NUCLEIC] [-O]
```

* `-p INPUT.pdb`
	: .pdb ファイル (入力)
* `-o OUTPUT.fred`
	: .fred ファイル (出力)
* `-sp FRAGMENT_OPTION_PROTEIN`, `--separate-protein FRAGMENT_OPTION_PROTEIN`
	: タンパク質のフラグメント分割方法 (`+amino`, `/amino`, `+peptide`, `/peptide`; デフォルト: `+amino`)
* `-sn FRAGMENT_OPTION_NUCLEIC`, `--separate-nucleic FRAGMENT_OPTION_NUCLEIC`
	: 核酸のフラグメント分割方法 (`+base`, `/base`, `/sugar`; デフォルト: `+base`)
* `-O`
	: 上書きプロンプトを表示しない


#### edit mode
```sh
$ fred4.py edit [-h] -i INPUT.ajf -o OUTPUT.fred [-p REF.pdb] [-O]
```

* `-i INPUT.ajf`
	: .ajf ファイル (入力)
* `-o OUTPUT.fred`
	: .fred ファイル (出力)
* `-p REF.pdb`
	: .pdb ファイル (指定無しの場合、ReadGeom のファイルを読み込む)
* `-O`
	: 上書きプロンプトを表示しない


#### rewrite mode
```sh
$ fred4.py rewrite [-h] -i INPUT.fred -o OUTPUT.fred [-O]
```

* `-i INPUT.fred`
	: .fred ファイル (入力)
* `-o OUTPUT.fred`
	: .fred ファイル (出力)
* `-O`
	: 上書きプロンプトを表示しない


#### output mode
```sh
$ fred4.py output [-h] -i INPUT.fred -o OUTPUT.ajf [-p REF.pdb] [-O]
```

* `-i INPUT`
	: .fred ファイル (入力)
* `-o OUTPUT`
	: .ajf ファイル (出力) (指定なしの場合、標準出力)
* `-p REF.pdb`
	: .pdb ファイル (指定無しの場合、ReadGeom を読み込む)
* `-O`
	: 上書きプロンプトを表示しない


#### writefrag mode
```sh
$  fred4.py writefrag [-h] -i INPUT.fred -p STRUCTURE.pdb -o OUTPUT_PREFIX
```

* `-i INPUT`
	: .fred ファイル (入力)
* `-p STRUCTURE.pdb`
	: 全体構造の .pdb ファイル
* `-o OUTPUT_PREFIX`
	: 出力ファイルの接頭辞


### editfrag.py
```sh
$ editfrag.py [-h] -b WHOLE.pdb -n ADD.pdb [ADD.pdb ...] -f BASE.fred -o NEW.fred [-c ATOM1-ATOM2 [ATOM1-ATOM2 ...]] [-m] [-O]
```

* `-b WHOLE.pdb`
	: 全体構造の .pdb ファイル
* `-n ADD.pdb`
	: フラグメントの .pdb ファイル
* `-f BASE.fred`
	: フラグメント分割前の .ajf ファイル
* `-o NEW.fred`
	: .fred ファイル (出力)
* `-c ATOM1-ATOM2 [ATOM1-ATOM2 ...]`
	: 接続情報を AMBERMASK で指定 (例: :LIG@C9-:LIG@C10)
* `-m`
	: 同じフラグメントの分割方法を他の同種の分子にも適用する。
* `-O`
	: 上書きプロンプトを表示しない




## fred ファイルの仕様
```txt
  FNo.  | Charge | BDA | Atoms of fragment
      1 |    1   |  0  |        1        2        3 …
      2 |   -1   |  1  |       11       12       13 …
      3 |    0   |  1  |       40       41       42 …
      4 |   -1   |  1  |       45       46       47 …
      5 |    0   |  1  |       74       75       76 …
      :      :      :

<< connections (ex. "Next_fragment_atom   Prev_fragment_atom") >>
        9        11
       37        40
       43        45
        :         :

===============< namelist >===============
&CNTRL
  Title='C9orf72-GA8_md-samp_11000_5A'
  ElecState='S1'
  Method='MP2'
  Nprint=3
        :
&FRAGMENT
{...}
/
```

* 1 行目: フラグメント情報の列名
* 2 行目〜 `<< connections >>` までの行
	* 1 行 1 フラグメントで記述
	* 各項目は `|` で識別している。
	* 1 列目
		* フラグメント番号
		* プログラム上、それほど重要ではなく、動作上はフラグメント情報の順番で識別している。
	* 2 列目
		* フラグメント電荷
	* 3 列目
		* BDA
		* 他のフラグメントへの接続数
		* 先頭フラグメントは 0。
	* 4 列目
		* フラグメントを構成する原子の原子インデックス
		* 固定長ではなく、連続スペースを区切り文字として識別している。
* `<< connections >>` 〜 `=====< namelist >=====` までの行
	* BAA
	* フラグメント間を接続している原子インデックス
	* 「次のフラグメントの原子インデックス」→「前のフラグメントの原子インデックス」の順に記述する。
* `=====< namelist >=====` 以降の行
	* .ajf ファイルの内容
	* この記述を変えると `output` モードで、変えた内容になる。
	* `&FRAGMENT` ネームリストは `{...}` で表現され、`output` モードではこの記号が置換対象となる。


## 動作要件
* python3
	* parmed


## License
The MIT License (MIT)

Copyright (c) 2021 Tatsuya Ohyama


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 13.0 (2024-11-19)
* `autofrag.py` を parmed を使った仕様に変更した。
* fred4.py に `autofrag.py` を `autofrag` モードとして統合し、`autofrag.py` を廃止した。
* `fragseparator.py` を parmed を使った仕様に変更した。
* fred4.py に `fragseparator.py` を `writefrag` モードとして統合し、`fragseparator.py` を廃止した。

### Ver. 12.0 (2024-10-08)
* `autofrag` サブコマンドを実装した。
* `autofrag.py` を廃止した。
* `FragmentData` クラスを parmed オブジェクトに対応した。

### Ver. 11.0 (2021-12-22)
* 公開した。
