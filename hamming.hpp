/**
 * @author Новиков А.В., nawww83@gmail.com.
 */

#pragma once

#include <iostream> // std::cout
#include <cassert>  // std::assert
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
        for (const auto& el : row) {
            std::cout << el << ", ";
        }
        std::cout << '\n';
    }
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
      {
         return { SymbolStatus::Erased, {} };
      }
      if( ( mStatus == SymbolStatus::Uninitialized ) || ( other.mStatus == SymbolStatus::Uninitialized ) )
      {
         return { SymbolStatus::Uninitialized, {} };
      }
      CodeElement result{ SymbolStatus::Normal, {} };
      for( int i = 0; i < N; ++i )
      {
         result.mSymbol[ i ] = mSymbol[ i ] ^ other.mSymbol[ i ];
      }
      return result;
   }
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
inline bool FormLeadBySwap( int i, Matrix< T >& H, int column_idx = -1, const std::vector< int >& columns = {} )
{
   int R = H.size();
   int N = H.at( 0 ).size();
   const int column = column_idx == -1 ? N - R + i : column_idx;
   auto main_element = H.at( i ).at( column );
   if( main_element != 0 )
      return true;
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
      return false;
   }
   for( int j = 0; j < R; ++j )
   {
      std::swap( H[ j ][ column ], H[ j ][ idx ] );
   }
   return true;
}

/**
 * Формирует систематическую проверочную матрицу по несистематической.
 */
template< typename T >
inline Matrix< T > MakeParityMatrixSystematic( const Matrix< T >& H, bool& is_ok,
                                               const std::vector< int >& columns = {} )
{
   is_ok = true;
   int R = H.size();
   int N = H.at( 0 ).size();
   auto result = H;
   // Form upper triangle.
   for( int i = R - 1; i >= 0; --i )
   {
      int idx = static_cast< int >( columns.size() ) < R ? -1 : columns.at( i );
      bool has_lead = FormLeadBySum( i, result, idx );
      if( !has_lead )
      {
         has_lead = FormLeadBySwap( i, result, idx, columns );
      }
      is_ok &= has_lead;
      for( int j = i - 1; j >= 0; --j )
      {
         int idx = static_cast< int >( columns.size() ) < R ? N + i - R : columns.at( i );
         const int reference = result.at( j ).at( idx );
         if( reference == 0 )
            continue;
         for( int k = 0; k < N; ++k )
            result[ j ][ k ] ^= result[ i ][ k ];
      }
   }
   // Form lower triangle.
   for( int i = 0; i < R; ++i )
   {
      for( int j = i + 1; j < R; ++j )
      {
         int idx = static_cast< int >( columns.size() ) < R ? N + i - R : columns.at( i );
         const int reference = result.at( j ).at( idx );
         if( reference == 0 )
            continue;
         for( int k = 0; k < N; ++k )
            result[ j ][ k ] ^= result[ i ][ k ];
      }
   }
   return result;
}

constexpr int power2( int x )
{
   if( x >= 0 )
      return int( 1 ) << x;
   else
      return 1;
}

/**
 * Расширенный код Хэмминга в систематической форме в режиме стирания ошибок.
 * Кодовое расстояние равно 4.
 * R - количество проверочных символов.
 */
template< int R, int M, typename T >
struct HammingExtended
{
   /**
    * Длина кода.
    */
   static constexpr int N = power2( R - 1 );
   
   /**
    * Количество информационных символов.
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
      {
         if( i == 0 )
         {
            mH.emplace_back( N, 1 );
         }
         else
         {
            mH.emplace_back( N, 0 );
         }
      }
      int deg = N / 2;
      for( int i = 1; i < R; ++i )
      {
         for( int j = 0; j < N; ++j )
         {
            mH[ i ][ j ] = ( ( ( j + 1 ) / deg ) % 2 ) == 1;
         }
         deg /= 2;
      }
      show_matrix(mH, "H in non-systematic form:");
      bool is_ok;
      mH = MakeParityMatrixSystematic( mH, is_ok );
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
            if( mH.at( i ).at( k ) == 0 )
               continue;
            element = element + a.at( k );
         }
         result.push_back( element );
      }
      return result;
   }

   /**
    * Вычислить синдром по принятому вектору (без стираний).
    */
   CodeWord< T, M > CalcSyndrome( const CodeWord< T, M >& v )
   {
      CodeWord< T, M > result;
      for( int i = 0; i < R; ++i )
      {
         CodeElement< T, M > element{ .mStatus = SymbolStatus::Normal, .mSymbol = {} };
         for( int k = 0; k < N; ++k )
         {
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
      // Определяем индексы стертых символов.
      std::vector< int > ids;
      for( int i = 0; i < N; ++i )
      {
         if( v.at( i ).mStatus == SymbolStatus::Erased )
            ids.push_back( i );
      }
      const int erased = ids.size();
      if( erased >= D )
      {
         return false;
      }
      // Выбираем часть проверочной матрицы - подматрицу.
      Matrix< int > selected;
      for( int j = 0; j < R; ++j )
      {
         selected.emplace_back( erased, 0 );
         for( int i = 0; auto idx : ids )
         {
            selected[ j ][ i++ ] = mH.at( j ).at( idx );
         }
      }
      // Упрощаем выбранную матрицу: максимально минимизируем количество единиц в каждой строке.
      // O(R^2) ... O(R^3)
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
      assert( good_rows.size() == erased );
      // Формируем строки из исходной матрицы в соответствие с индексами "хороших" строк.
      selected.clear();
      for( auto row : good_rows )
      {
         selected.push_back( mH.at( row ) );
      }
      // Делаем подматрицу систематической для прямого решения СЛАУ - восстановления стертых символов.
      {
         bool is_ok;
         selected = MakeParityMatrixSystematic( selected, is_ok, ids );
      }
      // Восстанавливаем стертые символы.
      for( int i = 0; i < erased; ++i )
      {
         CodeElement< T, M > recovered{ .mStatus = SymbolStatus::Normal, .mSymbol = {} };
         const int idx = ids.at( i );
         for( int k = 0; k < N; ++k )
         {
            if( k == idx )
               continue;
            if( selected.at( i ).at( k ) != 0 )
            {
               recovered = recovered + v.at( k );
            }
         }
         v[ idx ] = recovered;
      }
      return true;
   }

   /**
    * Проверочная матрица кода.
    */
   Matrix< int > mH;
};

template <typename T, int M>
inline void show_codeword(const CodeWord<T, M>& cword, const auto& code, const std::string& title) {
    std::cout << title << '\n';
    for (int k = 0; const auto& el1 : cword) {
        for (const auto& el2 : el1.mSymbol) {
            std::cout << int(el2) << ", ";
        }
        std::cout << '\n';
        if (k == (code.K - 1)) {
            std::cout << "----------" << '\n';
        }
        k++;
    }
    std::cout << '\n';
}

template <typename T, int M>
inline void show_cyndrome(const CodeWord<T, M>& c, const std::string& title) {
    std::cout << title << '\n';
    for (const auto& el1 : c) {
        for (const auto& el2 : el1.mSymbol) {
            std::cout << int(el2) << ", ";
        }
        std::cout << '\n';        
    }    
}

} // namespace hamming
