/**
 * @author Новиков А.В., nawww83@gmail.com.
 */

 #pragma once

 #include <iostream> // std::cout
 #include <cassert>  // assert
 #include <utility>  // std::pair
 #include <set>      // std::set
 #include <tuple>    // std::tie, std::ignore
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
     * Оператор сложения. В данном случае побитовый XOR.
     */
    CodeElement operator+( const CodeElement& other ) const
    {
       if( ( mStatus == SymbolStatus::Erased ) || ( other.mStatus == SymbolStatus::Erased ) )
          return { SymbolStatus::Erased, {} };
       if( ( mStatus == SymbolStatus::Uninitialized ) || ( other.mStatus == SymbolStatus::Uninitialized ) )
          return { SymbolStatus::Uninitialized, {} };
       CodeElement result{ SymbolStatus::Normal, {} };
       for( int i = 0; i < N; ++i )
          result.mSymbol[ i ] = mSymbol[ i ] ^ other.mSymbol[ i ];
       return result;
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
 inline std::pair<Matrix< T >, std::vector<std::pair<int, int>>> MakeParityMatrixSystematic( const Matrix< T >& H, bool& is_ok,
                                                const std::vector< int >& columns = {} )
 {
    is_ok = true; // Признак успешности преобразования.
    const int R = H.size();
    const int N = H.empty() ? 0 : H.at( 0 ).size();
    auto result = H;
    std::vector<std::pair<int, int>> swaps;
    if (H.empty()) {
      return std::make_pair(result, swaps);
    }
    std::pair<int, int> swaped_indexes;
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
 template< int R, int M, typename T >
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
     * Конструктор. Заполняется проверочная матрица кода, если не передана внешняя матрица.
     * H - Внешняя несистематическая проверочная матрица. Для задания произвольного кода.
     */
    explicit HammingExtended(const Matrix<int> H = {})
    {
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
         N = H.at(0).size();
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
    bool Decode( CodeWord< T, M >& v, int& erased )
    {
       assert(v.size() == N && "Input size is wrong");
       const auto& parity_check = mHsys;
       if (!mIsSystematic) {
          for (const auto& [a, b] : mSwaps) {
             std::swap( v[ a ], v[ b ] );
          }
       }
       // Определяем индексы стертых символов.
       std::vector< int > ids;
       for( int i = 0; i < N; ++i )
       {
          if( v.at( i ).mStatus == SymbolStatus::Erased )
             ids.push_back( i );
       }
       erased = ids.size();
       // Выбираем часть проверочной матрицы - подматрицу.
       Matrix< int > selected;
       for( int j = 0; j < R; ++j )
       {
          selected.emplace_back( erased, 0 );
          for( int i = 0; auto idx : ids )
             selected[ j ][ i++ ] = parity_check.at( j ).at( idx );
       }
      //  show_matrix(selected, "Selected columns:");
       // Упрощаем выбранную матрицу: максимально минимизируем количество единиц в каждой строке.
       // O(R^2).
       Matrix< int > selected_m = selected;
       selected = parity_check;
       Vector< int > row( erased, 0 );
       Vector< int > weights( R, -1 );
       // Метод 1. Дубовый. Раундовый.
      //  int counter = 0;
      //  for (;;) { // Раунды: итерации.
      //    counter++;
      //    bool all_minimal_weights = true;
      //    for( int i = 0; i < R; ++i )
      //    {
      //       int w = 0;
      //       for( int k = 0; k < erased; ++k )
      //       {
      //          w += selected_m.at( i ).at( k ) != 0; // Вес оригинального слова.
      //       }
      //       if (weights.at(i) == -1) {
      //          weights[i] = w;
      //       }
      //       for( int j = 0; j < R; ++j )
      //       {
      //          if( j == i )
      //             continue;
      //          int weight = 0;
      //          int bad_condition = 0;
      //          for( int k = 0; k < erased; ++k )
      //          {
      //             if ( selected_m.at( i ).at( k ) == 0 && selected_m.at( j ).at( k )) {
      //                bad_condition = 1; // Ноль в оригинале станет ненулевым элементом.
      //                break;
      //             }
      //             row[ k ] = selected_m.at( i ).at( k ) ^ selected_m.at( j ).at( k );
      //             weight += row[ k ] != 0;
      //          }
      //          if( weight < weights.at(i) && !bad_condition ) {
      //             selected_m[ i ] = row;
      //             weights[i] = weight;
      //             for( int k = 0; k < N; ++k )
      //             {
      //                selected[ i ][ k ] ^= selected.at( j ).at( k );
      //             }
      //          }
      //       }
      //    }
      //    for (auto w : weights) {
      //       all_minimal_weights &= (w <= 1);
      //    }
      //    if (all_minimal_weights || counter > erased) {
      //       break;
      //    }
      //  }
      //  show_matrix(selected_m, "Selected simplified:");
      //  std::cout << "counter: " << counter << std::endl;
       // Метод 2. Однопроходной (безраундовый).
       for( int i = 0; i < R; ++i )
       {
          for( int j = 0; j < R; ++j )
          {
             if( j == i )
                continue;
             int weight_original = 0;
             int weight = 0;
             for( int k = 0; k < erased; ++k )
             {
                row[ k ] = selected_m.at( i ).at( k ) ^ selected_m.at( j ).at( k );
                weight_original += selected_m.at( i ).at( k ) != 0;
                weight += row[ k ] != 0;
             }
             if( weight < weight_original ) {
                selected_m[ i ] = row;
                for( int k = 0; k < N; ++k )
               {
                  selected[ i ][ k ] ^= selected.at( j ).at( k );
               }
             }
          }
       }
      //  show_matrix(selected_m, "Selected simplified:");
       // Ищем базис - строки с единственной единицей ("хорошие" строки).
       std::vector< int > good_rows;
       std::set<std::vector<int>> uniques;
       for( int j = 0; j < R; ++j )
       {
         if (good_rows.size() == erased)
            break;
          int weight = 0;
          for( int k = 0; k < erased; ++k )
             weight += selected_m.at( j ).at( k ) != 0;
          if( weight == 1 && !uniques.contains(selected_m.at( j )) ) {
             good_rows.push_back( j );
             uniques.insert(selected_m.at( j ));
          }
       }
      //  show_matrix(selected, "Selected 1:");
       // Восстанавливаем стертые символы.
       for( const auto row : good_rows )
       {
         // std::cout << "i: " << i << ", M: " << M << ", v size: " << v.size() << std::endl;
          CodeElement< T, M > recovered{ .mStatus = SymbolStatus::Normal, .mSymbol = {} };
          int where_unit = -1;
            for (int k = 0; const auto el : selected_m.at( row )) {
               if (el) {
                  where_unit = k;
                  break;
               }
               k++;
            }
          for( int k = 0; k < N; ++k )
          {
             if( k == ids.at(where_unit) )
                continue;
             if( selected.at( row ).at( k ) != 0 ) {
               //  std::cout << "k: " << k << ", idx: " << idx << ", selected size: " << selected.size() << ", erased: " << erased << std::endl;
                assert(v.at(k).mStatus == SymbolStatus::Normal);
                recovered = recovered + v.at( k );
             }
          }
         //  std::cout << "idx: " << idx << ", erased: " << erased << ", good rows: " << good_rows.size() << std::endl;
          v[ ids.at(where_unit) ] = recovered;
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
    std::vector<std::pair<int, int>> mSwaps;
 
    /**
     * Проверочная матрица кода (несистематическая).
     */
    Matrix< int > mH;
 
    /**
     * Проверочная матрица систематического кода.
     */
    Matrix< int > mHsys;
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
 