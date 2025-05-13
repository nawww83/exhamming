/**
 * @author Новиков А.В., nawww83@gmail.com.
 */

#pragma once

#include <iostream> // std::cout
#include <cassert>  // assert
#include <utility>  // std::pair
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

template <typename T>
inline void show_vector(const Vector<T>& v, const std::string& title) {
   std::cout << title << '\n';
   for (const auto& el : v) {
      std::cout << el << ", ";
   }
   std::cout << std::endl;
}

template <typename T>
inline void show_vector(const Vector<std::pair<T, T>>& v, const std::string& title) {
   std::cout << title << '\n';
   for (const auto& [a, b] : v) {
      std::cout << "(" << a << ": " << b << "), ";
   }
   std::cout << std::endl;
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
   std::array< T, N > mSymbol;
   
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
inline bool FormLeadBySum( int i, Matrix< T >& H, std::vector< int >& correspondance, int column_idx = -1 )
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
   if (!correspondance.empty())
      correspondance[i] ^= correspondance.at(idx);
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
inline std::pair<Matrix< T >, std::vector<std::pair<int, int>>> MakeParityMatrixSystematic( const Matrix< T >& H, std::vector< int >& correspondance, bool& is_ok,
                                               const std::vector< int >& columns = {} )
{
   is_ok = true; // Признак успешности преобразования.
   const int R = H.size();
   const int N = H.empty() ? 0 : H.at( 0 ).size();
   auto result = H;
   std::vector<std::pair<int, int>> swaps;
   std::pair<int, int> swaped_indexes;
   if (H.empty()) {
      return std::make_pair(result, swaps);
   }
   // Формирование верхней треугольной матрицы (справа).
   for( int i = R - 1; i >= 0; --i )
   {
      const int idx = std::cmp_not_equal(columns.size(), R) ? -1 : columns.at( i );
      bool has_lead = FormLeadBySum( i, result, correspondance, idx );
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
         if (!correspondance.empty())
            correspondance[j] ^= correspondance.at(i);
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
         if (!correspondance.empty())
            correspondance[j] ^= correspondance.at(i);
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
   static constexpr auto N = power2<int>( R - 1 );
   
   /**
    * Количество информационных кодовых символов.
    */
   static constexpr int K = N - R;
   
   /**
    * Кодовое расстояние.
    */
   static constexpr int D = 4;

   /**
    * Конструктор. Заполняется проверочная матрица кода.
    */
   explicit HammingExtended()
   {
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
      bool is_ok;
      Vector<int> correspondance;
      std::tie(mHsys, mSwaps) = MakeParityMatrixSystematic( mH, correspondance, is_ok );
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
      auto& parity_check = mIsSystematic ? mHsys : mH;
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
   bool Decode( CodeWord< T, M >& v )
   {
      assert(v.size() == N);
      auto& parity_check = mIsSystematic ? mHsys : mH;
      // Определяем индексы стертых символов.
      std::vector< int > ids;
      for( int i = 0; i < N; ++i )
      {
         if( v.at( i ).mStatus == SymbolStatus::Erased )
            ids.push_back( i );
      }
      const int erased = ids.size();
      // Выбираем часть проверочной матрицы - подматрицу.
      Matrix< int > selected;
      for( int j = 0; j < R; ++j )
      {
         selected.emplace_back( erased, 0 );
         for( int i = 0; auto idx : ids )
            selected[ j ][ i++ ] = parity_check.at( j ).at( idx );
      }
      // Упрощаем выбранную матрицу: максимально минимизируем количество единиц в каждой строке.
      // O(R^2).
      Matrix< int > selected_m = selected;
      Vector< int > row( erased, 0 );
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
            if( weight < weight_original )
            {
               selected_m[ i ] = row;
            }
         }
      }
      // Ищем базис - строки с единственной единицей ("хорошие" строки).
      std::vector< int > good_rows;
      for( int j = 0; j < R; ++j )
      {
         int weight = 0;
         for( int k = 0; k < erased; ++k )
         {
            weight += selected_m.at( j ).at( k ) != 0;
         }
         if( weight == 1 )
         {
            good_rows.push_back( j );
         }
      }
      // assert( good_rows.size() == erased );
      // Формируем строки из исходной матрицы в соответствие с индексами "хороших" строк.
      selected.clear();
      for( auto row : good_rows )
         selected.push_back( parity_check.at( row ) );
      // Делаем подматрицу систематической для прямого решения СЛАУ - восстановления стертых символов.
      bool is_ok;
      Vector<int> correspondance;
      std::tie(selected, std::ignore) = MakeParityMatrixSystematic( selected, correspondance, is_ok, ids );
      // assert( selected.size() == erased );
      // Восстанавливаем стертые символы.
      for( int i = 0; i < std::min(selected.size(), ids.size()); ++i )
      {
         CodeElement< T, M > recovered{ .mStatus = SymbolStatus::Normal, .mSymbol = {} };
         const int idx = ids.at( i );
         for( int k = 0; k < N; ++k )
         {
            if( k == idx )
               continue;
            if( selected.at( i ).at( k ) != 0 )
               recovered = recovered + v.at( k );
         }
         v[ idx ] = recovered;
      }
      if (!mIsSystematic) {
         for (const auto& [a, b] : mSwaps) {
            std::swap( v[ a ], v[ b ] );
         }
      }
      while (v.size() > K)
         v.pop_back();
      return true;
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
inline void show_codeword(const CodeWord<T, M>& cword, const auto& code, const std::string& title) {
   std::cout << title << '\n';
   for (int k = 0; const auto& el1 : cword) {
      for (const auto& el2 : el1.mSymbol)
         std::cout << int(el2) << ", ";
      std::cout << '\n';
      if (k == (code.K - 1))
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
