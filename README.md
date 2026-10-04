# cppcplib

Nim版 `cplib/src/cplib` 全306ソースのC++20移植です。

長さ・要素数を持つ型では `.size()` を利用できます。既存の `.len()` と自由関数
`len(s)` はそのまま利用できます。返り値は `cplib::Int`、const オブジェクトにも対応します。

```cpp
#include <cplib/str/static_string.hpp>

auto s = cplib::toStaticString("abc");
auto n = s.size(); // 3
auto part = s.substr(1, 2);
auto m = part.size(); // 1
```

Graph の `.size()` は頂点数、UnionFind 類は構築時の総要素数です。辺数や連結成分数の
既存APIは変わりません。BitSet / BitVector は論理ビット数、WordsizeTree 類は登録要素数を
返します。オンライン畳み込みは確定済み係数数、遡及優先度キューは操作後に残る要素数です。
結合文字列は既存 `.len()` と同じ計算量です（結合数に比例する版もあります）。

HashString、RepeatedStaticString、StaticStringBase の公開 `size` フィールドは互換性のため
維持しています。これらの型には同名の `.size()` を追加していません。行列・二次元木の
次元や RangeSet の区間数/区間内要素数、回文木の文字数/節点数は一律に変更していません。

`sheep.hpp` の `len(x)` マクロは `s.len()` にも展開されます。ライブラリヘッダーを先に
読み込み、テンプレートを後に読み込む構成では `.size()` または `len(s)` を使えます。
メンバ/自由関数を直接呼ぶ場合は `(s.len)()`、`(cplib::len)(s)` でマクロ展開を避けられます。
