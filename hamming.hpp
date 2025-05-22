/**
 * @author Новиков А.В., nawww83@gmail.com.
 */

 #pragma once

 #include <iostream> // std::cout
 #include <cassert>  // assert
 #include <utility>  // std::pair
 #include <tuple>    // std::tie
 #include <array>    // std::array
 #include <vector>   // std::vector
 #include <string>   // std::string
 
 namespace hamming
 {
 
 template< typename T >
 using Vector = std::vector< T >;
 
 template< typename T >
 using Matrix = Vector< Vector< T > >;

 template <typename T>
 using Swaps = Vector<std::pair<T, T>>;
 
 template <typename T>
 inline void show_matrix(const Matrix<T>& M, const std::string& title) {
    std::cout << title << '\n';
    for (const auto& row : M) {
       for (const auto& el : row)
          std::cout << el << ", ";
       std::cout << '\n';
    }
    std::cout << std::flush;
 }
 
 /**
  * Статус принятого (канального) символа.
  */
 enum class SymbolStatus
 {
    Uninitialized = 0, // Неопределен.
    Normal,            // Обычное состояние.
    Erased             // Стертый.
 };
 
 /**
  * Кодовый элемент (символ).
  */
 template< typename T, int N >
 struct CodeElement
 {
    /**
     * Статус кодового символа.
     */
    SymbolStatus mStatus;
 
    /**
     * Внутренние символы.
     */
    std::array< T, N > mSymbol{};
    
    /**
     * Оператор сложения. В данном случае побитовый XOR. Умно комбинирует стертый символ с нестертым.
     */
    CodeElement operator+( const CodeElement& other ) const
    {
       if( ( mStatus == SymbolStatus::Normal ) && ( other.mStatus == SymbolStatus::Normal ) )
       {
         CodeElement result{ SymbolStatus::Normal, {} };
         for( int i = 0; i < N; ++i )
            result.mSymbol[ i ] = mSymbol[ i ] ^ other.mSymbol[ i ];
         return result;
       }
       if( ( mStatus == SymbolStatus::Erased ) && ( other.mStatus == SymbolStatus::Normal ) )
          return { SymbolStatus::Normal, other.mSymbol }; // 0 + x = 0; 1 + x = 1.
       if( ( mStatus == SymbolStatus::Normal ) && ( other.mStatus == SymbolStatus::Erased ) )
          return { SymbolStatus::Normal, mSymbol };
       if( ( mStatus == SymbolStatus::Erased ) && ( other.mStatus == SymbolStatus::Erased ) )
          return { SymbolStatus::Erased, {} }; // x + x = x.
      return { SymbolStatus::Uninitialized, {} }; 
    }
 
    /**
     * 
     */
    bool operator==( const CodeElement< T, N >& other ) const = default;
 };
 
 /**
  * Кодовое слово (вектор).
  */
 template< typename T, int N >
 using CodeWord = std::vector< CodeElement< T, N > >;
 
 /**
  * Формирует лидирующие элементы, используя взвешенную сумму.
  * Таким образом у проверочной матрицы справа формируется единичная матрица.
  */
 template< typename T >
 inline bool FormLeadBySum( int i, Matrix< T >& H, int column_idx = -1 )
 {
    assert(!H.empty());
    int R = H.size();
    int N = H.at( 0 ).size();
    const int column = column_idx == -1 ? N - R + i : column_idx;
    auto main_element = H.at( i ).at( column );
    if( main_element != 0 )
       return true;
    int idx = -1;
    for( int j = i - 1; j >= 0; --j )
    {
       auto element = H.at( j ).at( column );
       if( element != 0 )
       {
          idx = j;
          break;
       }
    }
    if( idx == -1 )
    {
       return false;
    }
    for( int k = 0; k < N; ++k )
    {
       H[ i ][ k ] ^= H.at( idx ).at( k );
    }
    return true;
 }
 
 /**
  * Формирует лидирующие элементы, используя перестановки столбцов (swap).
  * Таким образом у проверочной матрицы справа формируется единичная матрица.
  */
 template< typename T >
 inline std::pair<bool, std::pair<int, int>> FormLeadBySwap( int i, Matrix< T >& H, int column_idx = -1, const std::vector< int >& columns = {} )
 {
   assert(!H.empty());
    int R = H.size();
    int N = H.at( 0 ).size();
    const int column = column_idx == -1 ? N - R + i : column_idx;
    auto main_element = H.at( i ).at( column );
    if( main_element != 0 )
       return std::make_pair(true, std::make_pair(-1, -1));
    int idx = -1;
    if( columns.size() == 0 )
    {
       for( int j = 0; j < N - R; ++j )
       {
          auto element = H.at( i ).at( j );
          if( element != 0 )
          {
             idx = j;
             break;
          }
       }
    }
    else
    {
       for( int j = 0; j < N; ++j )
       {
          bool has = false;
          for( auto el : columns )
          {
             has |= el == j;
          }
          if( has )
             continue;
          auto element = H.at( i ).at( j );
          if( element != 0 )
          {
             idx = j;
             break;
          }
       }
    }
    if( idx == -1 )
    {
       return std::make_pair(false, std::make_pair(-1, -1));
    }
    for( int j = 0; j < R; ++j )
    {
       std::swap( H[ j ][ column ], H[ j ][ idx ] );
    }
    return std::make_pair(true, std::make_pair(column, idx));
 }

 /**
  * Формирует систематическую проверочную матрицу по несистематической.
  * @param columns Столбцы, которые будут базисными (по умочанию - справа).
  */
 template< typename T >
 inline std::pair<Matrix< T >, Swaps<T>> MakeParityMatrixSystematic( const Matrix< T >& H, bool& is_ok,
                                                const std::vector< int >& columns = {} )
 {
    is_ok = true; // Признак успешности преобразования.
    const int R = H.size();
    const int N = H.empty() ? 0 : H.at( 0 ).size();
    auto result = H;
    Swaps<T> swaps;
    if (H.empty()) {
      return std::make_pair(result, swaps);
    }
    std::pair<T, T> swaped_indexes;
    // Формирование верхней треугольной матрицы (справа).
    for( int i = R - 1; i >= 0; --i )
    {
       const int idx = std::cmp_not_equal(columns.size(), R) ? -1 : columns.at( i );
       bool has_lead = FormLeadBySum( i, result, idx );
       if( !has_lead ) {
          std::tie(has_lead, swaped_indexes) = FormLeadBySwap( i, result, idx, columns );
          swaps.push_back(swaped_indexes);
       }
       is_ok &= has_lead;
       for( int j = i - 1; j >= 0; --j )
       {
          int idx = std::cmp_not_equal(columns.size(), R) ? N + i - R : columns.at( i );
          if( auto reference = result.at( j ).at( idx ); reference == 0 )
             continue;
          for( int k = 0; k < N; ++k )
             result[ j ][ k ] ^= result[ i ][ k ];
       }
    }
    // Формирование нижней треугольной матрицы (справа).
    for( int i = 0; i < R; ++i )
    {
       for( int j = i + 1; j < R; ++j )
       {
          int idx = std::cmp_not_equal(columns.size(), R) ? N + i - R : columns.at( i );
          if( auto reference = result.at( j ).at( idx ); reference == 0 )
             continue;
          for( int k = 0; k < N; ++k )
             result[ j ][ k ] ^= result[ i ][ k ];
       }
    }
    return std::make_pair(result, swaps);
 }

 /**
  * Упрощает проверочную матрицу parity_check, чтобы решить СЛАУ относительно неизвестных, 
  * которые соответствуют столбцам selected из проверочной матрицы.
  * O(R^2) * E, где E - количество стертых символов.
  */
template< typename T, int M >
void SimplifyManyOfPass(std::vector< CodeElement<T, M> >& free_column, Matrix<int>& selected, int swap_direction) {
   [[maybe_unused]] int counter = 0;
   const int R = free_column.size();
   assert(R > 0);
   assert(!selected.empty());
   const int erased = selected.at(0).size();
   assert(erased > 0);
   Vector< int > row( erased );
   Vector< int > weights( R);
   // Проход "сверху вниз".
   auto up_down = [&]() -> bool {
      bool no_modifications = true;
      for( int i = 0; i < R; ++i )
      {
         int w = 0;
         for( int k = 0; k < erased; ++k )
            w += selected.at( i ).at( k ) != 0; // Вес оригинального слова.
         weights[i] = w;
         for( int j = i+1; j < R; ++j )
         {
            int weight = 0;
            for( int k = 0; k < erased; ++k )
            {
               row[ k ] = selected.at( i ).at( k ) ^ selected.at( j ).at( k );
               weight += row[ k ] != 0;
            }
            if( weight < weights.at(i) ) {
               no_modifications = false;
               selected[ i ] = row;
               weights[i] = weight;
               free_column[ i ] = free_column.at(i) + free_column.at( j );
            }
         }
      }
      return no_modifications;
   };
   // Проход "снизу вверх".
   auto down_up = [&]() -> bool {
      bool no_modifications = true;
      for( int i = R-1; i >= 0; --i )
      {
         int w = 0;
         for( int k = 0; k < erased; ++k )
            w += selected.at( i ).at( k ) != 0; // Вес оригинального слова.
         weights[i] = w;
         for( int j = i-1; j >= 0; --j )
         {
            int weight = 0;
            for( int k = 0; k < erased; ++k )
            {
               row[ k ] = selected.at( i ).at( k ) ^ selected.at( j ).at( k );
               weight += row[ k ] != 0;
            }
            if( weight < weights.at(i) ) {
               no_modifications = false;
               selected[ i ] = row;
               weights[i] = weight;
               free_column[ i ] = free_column.at(i) + free_column.at( j );
            }
         }
      }
      return no_modifications;
   };
   for (;;) { // Раунды: итерации.
      counter++;
      bool no_modifications_1;
      bool no_modifications_2;
      // std::cout << "Counter " << counter << std::endl;
      if (swap_direction == 0) {
         no_modifications_1 = up_down();
         no_modifications_2 = down_up();
      } else {
         no_modifications_2 = down_up();
         no_modifications_1 = up_down();
      }
      if (no_modifications_1 && no_modifications_2)
         break;
   }
}

 template <typename T>
 inline constexpr T power2( int x )
 {
    return (x > 0) ? (T( 1 ) << x) : T(1);
 }
 
 /**
  * Расширенный векторный код Хэмминга. Декодирование в режиме стирания ошибок.
  * R - количество проверочных символов.
  * M - количество внутренних символов в одном кодовом символе (векторность кода).
  * T - тип внутреннего символа. Для упрощения может быть шире типа реально используемых данных.
  */
 template<typename T, int R, int M>
 struct HammingExtended
 {
    /**
     * Длина кода. Код двоичный в плане кодового символа, однако, кодовый символ - векторный.
     * Используется алгебра побитового XOR, которой нет разницы сколько внутренних символов (элементов вектора).
     */
    int N = power2<int>( R - 1 );
    
    /**
     * Количество информационных кодовых символов.
     */
    int K = N - R;

    /**
     * Кодовое расстояние кода. По умолчанию - расстояние для расширенного кода Хэмминга.
     */
    int D = 4;
 
    /**
     * Конструктор. Заполняется проверочная матрица кода, если не передана внешняя матрица.
     * H - Внешняя несистематическая проверочная матрица. Для задания произвольного кода.
     */
    explicit HammingExtended(const Matrix<int> H = {}, int code_distance = -1)
    {
      assert((!H.empty() && code_distance != -1) || (H.empty() && code_distance == -1));
      if (code_distance != -1) {
         D = code_distance;
      }
      assert(D <= (R + 1));
      if (H.empty()) {
         mH.clear();
         for( int i = 0; i < R; ++i )
            mH.emplace_back( N, i == 0 );
         int deg = N / 2;
         for( int i = 1; i < R; ++i )
         {
            for( int j = 0; j < N; ++j )
               mH[ i ][ j ] = ( ( ( j + 1 ) / deg ) % 2 ) == 1;
            deg /= 2;
         }
      } else {
         assert(H.size() == R);
         assert(R > 0);
         N = H.at(0).size();
         assert(N > R);
         K = N - R;
         mH = H;
      }
      bool is_ok;
      std::tie(mHsys, mSwaps) = MakeParityMatrixSystematic( mH, is_ok );
      // show_matrix(mHsys, "Systematic:");
      assert(is_ok);
    }
 
    /**
     * Закодировать информационный вектор.
     */
    CodeWord< T, M > Encode( const CodeWord< T, M >& a )
    {
       CodeWord< T, M > result;
       for( const auto& el : a )
       {
          assert( el.mStatus == SymbolStatus::Normal );
          result.push_back( el );
       }
       for( int i = 0; i < R; ++i )
       {
          CodeElement< T, M > element{ .mStatus = SymbolStatus::Normal, .mSymbol = {} };
          for( int k = 0; k < K; ++k )
          {
             if( mHsys.at( i ).at( k ) == 0 )
                continue;
             element = element + a.at( k );
          }
          result.push_back( element );
       }
       if (!mIsSystematic) {
          for (const auto& [a, b] : mSwaps) {
             std::swap( result[ a ], result[ b ] );
          }
       }
       return result;
    }
 
    /**
     * Вычислить синдром по принятому вектору (без стираний).
     */
    CodeWord< T, M > CalcSyndrome( const CodeWord< T, M >& v )
    {
       CodeWord< T, M > result;
       const auto& parity_check = mIsSystematic ? mHsys : mH;
       for( int i = 0; i < R; ++i )
       {
          CodeElement< T, M > element{ .mStatus = SymbolStatus::Normal, .mSymbol = {} };
          for( int k = 0; k < N; ++k ) {
             if (parity_check.at(i).at(k) == 0)
                continue;
             element = element + v.at( k );
          }
          result.push_back( element );
       }
       return result;
    }
 
    /**
     * Декодировать принятый вектор в режиме стирания ошибки.
     */
    bool Decode( CodeWord< T, M >& v, int& erased, int& was_changed_strategy )
    {
       assert(v.size() == N && "Input size is wrong");
       if (!mIsSystematic) {
          for (const auto& [a, b] : mSwaps) {
             std::swap( v[ a ], v[ b ] );
          }
       }
       const auto& parity_check = mHsys;
       // Определяем индексы стертых символов и формируем столбец свободных членов.
       std::vector< int > ids;
       mFreeColumn.clear();
       for (int j = 0; j < R; ++j) {
         mFreeColumn.push_back(CodeElement<T, M>{.mStatus = hamming::SymbolStatus::Normal, .mSymbol = {}});
       }
       for( int i = 0; i < N; ++i )
       {
         if( v.at( i ).mStatus == SymbolStatus::Erased ) {
            ids.push_back( i );
         } else {
            for( int j = 0; j < R; ++j )
            {
               if (parity_check.at(j).at(i))
                  mFreeColumn[j] = mFreeColumn.at(j) + v.at(i);
            }
         }
       }
       erased = ids.size();
       if (erased > R) {
         return false;
       }
       // Выбираем часть проверочной матрицы - подматрицу.
       auto select_erasure_submatrix = [this, &parity_check](std::vector< int >& ids) {
          mErasureSubmatrix.clear();
          const int erased = ids.size();
         for( int j = 0; j < R; ++j )
         {
            mErasureSubmatrix.emplace_back( erased, 0 );
            for( int i = 0; auto idx : ids )
               mErasureSubmatrix[ j ][ i++ ] = parity_check.at( j ).at( idx );
         }
       };
       select_erasure_submatrix(ids);
       // Упрощаем подматрицу: минимизируем количество единиц в каждой строке.
       // Параллельно модифицируем столбец свободных членов.
       auto free_column = mFreeColumn;
       SimplifyManyOfPass(free_column, mErasureSubmatrix, 0);
       // Ищем базис по стертым столбцам - строки с единственной единицей - "хорошие" строки.
       std::vector< int > good_rows;
       auto select_good_rows = [this, &good_rows](int erased) -> int {
         int bad_strategy = 0;
         good_rows.clear();
         for( int j = 0; j < R; ++j )
         {
            if (good_rows.size() == erased)
               break;
            int weight = 0;
            for( int k = 0; k < erased; ++k ) {
               weight += mErasureSubmatrix.at( j ).at( k ) != 0;
            }
            if( weight == 1 ) {
               good_rows.push_back( j );
            }
            if (weight > 1) {
               bad_strategy = 1; // Выбрана неудачная стратегия упрощения матрицы.
            }
         }
         return bad_strategy;
       };
       int bad_strategy = select_good_rows(erased);
       was_changed_strategy = 0;
       if (bad_strategy) { // Обычно для больших стираний: около D-1 и выше.
         free_column = mFreeColumn;
         select_erasure_submatrix(ids);
         was_changed_strategy = 1;
         SimplifyManyOfPass(free_column, mErasureSubmatrix, 1); // Пробуем инвертированную стратегию.
         bad_strategy = select_good_rows(erased);
         // std::cout << "Change direct. Strategy after is " << (bad_strategy ? "bad" : "good") << ", erased: " << erased << std::endl;
       }
       // Восстанавливаем стертые символы.
       for( const auto& row : good_rows )
       {
         if (bad_strategy) { // Не нашли стратегию: считаем, что система однозначно не разрешима.
            break;
         }
         // std::cout << "i: " << i << ", M: " << M << ", v size: " << v.size() << std::endl;
          CodeElement< T, M > recovered{ .mStatus = SymbolStatus::Normal, .mSymbol = {} };
          int where_unit = -1;
         for (int k = 0; const auto& el : mErasureSubmatrix.at( row )) {
            if (el) {
               where_unit = k;
               break;
            }
            k++;
         }
         //  std::cout << "idx: " << idx << ", erased: " << erased << ", good rows: " << good_rows.size() << std::endl;
          v[ ids.at(where_unit) ] = free_column.at(row);
       }
       while (v.size() > K)
          v.pop_back();
       return true;
    }

    /**
     * 
     */
    auto getSwaps() const {
      return mSwaps;
    }
 
    /**
     * Включить/выключить режим систематического кодирования. 
     */
    void SwitchToSystematic(bool is_systematic) {
       mIsSystematic = is_systematic;
    }
 
    bool mIsSystematic = true;
 
    /**
     * Сделанные во время формирования систематической матрицы перестановки столбцов. 
     */
    Swaps< int > mSwaps;
 
    /**
     * Проверочная матрица кода (несистематическая).
     */
    Matrix< int > mH;
 
    /**
     * Проверочная матрица систематического кода.
     */
    Matrix< int > mHsys;

    /**
     * Выборочная матрица, соответствующая столбцам со стираниями.
     */
    Matrix< int > mErasureSubmatrix;

    /**
     * Столбец свободных членов.
     */
    std::vector<CodeElement<T, M>> mFreeColumn;
 };
 
 template <typename T, int M>
 inline void show_codeword(const CodeWord<T, M>& cword, int K, const std::string& title) {
    std::cout << title << '\n';
    for (int k = 0; const auto& el1 : cword) {
       for (const auto& el2 : el1.mSymbol)
          std::cout << int(el2) << ", ";
       std::cout << '\n';
       if (k == (K - 1))
          std::cout << "----------\n";
       k++;
    }
    std::cout << std::endl;
 }
 
 template <typename T, int M>
 inline void show_cyndrome(const CodeWord<T, M>& c, const std::string& title) {
    std::cout << title << '\n';
    for (const auto& el1 : c) {
       for (const auto& el2 : el1.mSymbol)
          std::cout << int(el2) << ", ";
       std::cout << '\n';        
    }
    std::cout << std::flush;
 }
 
 } // namespace hamming
 
