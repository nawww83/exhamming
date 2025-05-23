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
  * Подготавливает методом Гаусса матрицу и столбец свободных членов к решению обратным ходом.
  * @param free_column - столбец свободных членов, dim(free_column) = (r x 1).
  * @param selected - прямоугольная матрица, dim(selected) = (r x e).
  * Сложность O(r * e^2).
  */
 template< typename T, int M >
 void Gauss(CodeWord<T, M>& free_column, Matrix<int>& selected) {
   const int R = free_column.size();
   assert(R > 0);
   assert(!selected.empty());
   const int erased = selected.at(0).size();
   assert(erased > 0);
   for( int k = 0; k < erased; ++k ) {
      int where_unit = -1;
      for( int i = k; i < R; ++i ) {
         const bool is_not_zero = (selected.at(i).at(k) != 0);
         if (is_not_zero) {
            where_unit = i;
            break;
         }
      }
      if (where_unit == -1)
         continue;
      if (where_unit > k) {
         free_column[ k ] = free_column.at(k) + free_column.at( where_unit );
         for( int k1 = 0; k1 < erased; ++k1 ) {
            selected[k][k1] ^= selected.at( where_unit ).at( k1 );
         }
      }
      // Обнуляем до конца: при этом метод гарантированно за один проход выдает разрешимую 
      // матрицу (если исходная СЛАУ имеет решение).
      // Матрица разрешимая - значит все элементы главной диагонали квадратной подматрицы ненулевые.
      // При такой стратегии все элементы ниже квадратной подматрицы равны нулю.
      for( int i = k + 1; i < R; ++i ) {
         const bool is_not_zero = (selected.at(i).at(k) != 0);
         if (is_not_zero) {
            free_column[ i ] = free_column.at(i) + free_column.at( k );
            for( int k1 = 0; k1 < erased; ++k1 ) {
               selected[i][k1] ^= selected.at( k ).at( k1 );
            }
         }
      }
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
      auto free_column = mFreeColumn;
      // show_codeword(free_column, -1, "free column:");
      // show_matrix(mErasureSubmatrix, "Erasure matrix:");
      Gauss(free_column, mErasureSubmatrix);
      // show_codeword(free_column, -1, "free column after Gauss method:");
      // show_matrix(mErasureSubmatrix, "Erasure matrix after Gauss method:");
      // Восстанавливаем стертые символы: решение СЛАУ обратным ходом.
      for( int k = erased - 1; k >= 0; --k ) {
         // std::cout << " k = " << k << std::endl;
         const int idx_v = ids.at(k);
         if (mErasureSubmatrix.at(k).at(k) != 0) {
            v[ idx_v ] = free_column.at(k);
            // std::cout << " v = " << std::endl;
         }
         for (int j = 0; j < erased - 1 - k; j++) {
            if (mErasureSubmatrix.at(k).at(k + j + 1) != 0) {
               v[idx_v] = v.at(idx_v) + v.at(ids.at(k + j + 1));
               // std::cout << " v += " << std::endl;
            }
         }
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
    CodeWord<T, M> mFreeColumn;
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
 