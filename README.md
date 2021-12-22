# FragmentTools

## 概要
ABINIT-MP による FMO 計算で、手動フラグメント分割をサポートするツール群


## ツール群
* `fred4.py`
	: .ajf → .fred 間、.fred → .fred 間、.fred → .ajf 間でファイルを相互変換するプログラム。
	: `edit`、`rewrite`、`output`、`autofrag`、`editfrag` の 5 つの機能を持つ。
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
* `autofrag.py`
	: PDB ファイルから、自動でフラグメントを分割するプログラム。
* `editfrag.py`
	: フラグメント情報を構造を与えてフラグメントを編集する。
	: 複数の分子に対して同じフラグメント分割をする。


## 使用方法
### autofrag.py
```sh
$ autofrag.py [-h] -i INPUT.pdb -o OUTPUT.fred [-S] -V VERSION [-O]
```

* `-i INPUT.pdb`:
	: フラグメント分割をする PDB ファイル
* `-o OUTPUT.fred`
	: 人間可読な編集ファイルの出力先
* `-S`
	: DNA をリン酸基と分割するかどうか (デフォルト: OFF)
* `-V VERSION`:
	: .ajf ファイルのバージョンか、テンプレートファイル
	* `3`: ABINIT-MP 3
	* `5`: ABINIT-MP 5 or later
	* `m`: mizuho ABINIT-MP
* `-O`
	: 上書きプロンプトを表示しない


### fred4.py
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


#### output
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
### Ver. 11.0 (2021-12-22)
* 公開した。
