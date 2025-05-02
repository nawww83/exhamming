#include <iostream> // std::cout
#include <random>
#include <set>
#include <chrono>
#include "hamming.hpp"

namespace {
namespace rng_n {
   auto G_SEED = std::random_device{}();
   
   template <typename T>
   class GeometricDistribution {
   private:
      std::mt19937 gen;
      std::geometric_distribution<T> dist;
   public:
      GeometricDistribution(double p): gen(G_SEED), dist(p) {}
      GeometricDistribution(GeometricDistribution&& other) noexcept: dist(std::move(other.dist)), gen(std::move(other.gen)) {}
      T operator()() { return dist(gen); }
      void seed() {
         auto t = std::chrono::high_resolution_clock::now().time_since_epoch().count();
         std::seed_seq seq {t & 255, (t >> 8) & 255, (t >> 16) & 255, (t >> 24) & 255};
         gen.seed( seq );
      }
   };

   template <typename T>
   class UniformIntDistribution {
   private:
      std::mt19937 gen;
      std::uniform_int_distribution<T> dist;
   public:
      UniformIntDistribution(): gen(G_SEED) {}
      UniformIntDistribution(UniformIntDistribution&& other) noexcept: dist(std::move(other.dist)), gen(std::move(other.gen)) {}
      T operator()() { return dist(gen); }
      void seed() {
         auto t = std::chrono::high_resolution_clock::now().time_since_epoch().count();
         std::seed_seq seq {t & 255, (t >> 8) & 255, (t >> 16) & 255, (t >> 24) & 255};
         gen.seed( seq );
      }
   };
   
   template<int n, int m, class Generator, class Generator2>
   auto get_random_row(Generator& g, Generator2& g2) {
      std::vector<int> row(n, 0);
      int sum = 0;
      int idx = std::abs(g2()) % n;
      for (;;) {
         if (sum == m) {
            break;
         }
         for (;;) {
            if ((g() % 2) && !row.at(idx)) {
               row[idx] = 1;
               sum++;
            }
            idx++;
            idx %= n;
            if (sum == m) {
               break;
            }
         }
      }
      return row;
   }
}
}

static void TestHamming()
{
   using namespace hamming;
   using u8 = uint8_t;  // Тип внутреннего символа. В данном случае - байт.
   constexpr int R = 5; // Количество проверочных кодовых символов.
   constexpr int M = 4; // Количество внутренних символов на один кодовый символ (векторность).
   HammingExtended< R, M, u8 > code;
   code.SwitchToSystematic(true); // Cистематический код.
   // code.SwitchToSystematic(false); // Несистематический код.
   std::cout << "K: " << code.K << ", N: " << code.N << '\n';    
   show_matrix(code.mH, "Parity check matrix H:");
   show_matrix(code.mHsys, "Parity check matrix in systematic form H:");
   // Формирование некоторого информационного вектора.
   CodeWord< u8, M > a( code.K, { .mStatus = SymbolStatus::Normal, .mSymbol = {} } );
   for( u8 i = 0; auto& el : a )
   {
      std::array< u8, M > values{ i++, i++, i++, i++ };
      el.mSymbol = values;
   }
   // Кодирование.
   auto s = code.Encode( a );
   show_codeword(s, code, "Codeword:"); 
   // Вычисление синдрома.
   const auto& c = code.CalcSyndrome( s );
   show_cyndrome(c, "Cyndrome of the codeword:"); 
   // Стирание нескольких символов в канале. Кодовое расстояние кода равно 4,
   // кратность стирания q <= d - 1 = 3. Значит все стирания из не более трех 
   // символов будут успешно восстановлены.
   auto v = s;
   v[ 5 ] = { .mStatus = SymbolStatus::Erased, .mSymbol = {} };
   v[ 7 ] = { .mStatus = SymbolStatus::Erased, .mSymbol = {} };
   v[ 12 ] = { .mStatus = SymbolStatus::Erased, .mSymbol = {} };
   // v[ 2 ] = { .mStatus = SymbolStatus::Erased, .mSymbol = {} };
   // Декодирование (восстановление стертых символов).
   const auto decode_is_ok = code.Decode( v );
   const bool recover_is_ok = v == a;
   show_codeword(v, code, "Decoded symbols:");  
   std::cout << "Decode is : " << (decode_is_ok ? "Ok\n" : "Failed\n");
   std::cout << "Recover is : " << (recover_is_ok ? "Ok\n" : "Failed\n");
   std::cout << '\n';
}


// void MakeLdpc() {
//    hamming::Matrix<int> H;
//    hamming::Matrix<int> Hsys;
//    hamming::Vector<int> cyndrome;
//    std::vector<std::pair<int, int>> swaps;
//    rng_n::GeometricDistribution<int> gd(0.3);
//    gd.seed();
//    rng_n::UniformIntDistribution<int> ud;
//    ud.seed();
//    /**
//     * Длина LDPC кода. Для примера рассматривается относительно короткий код.
//     */
//    constexpr int n = 32;
//    /**
//     * Количество проверочных символов.
//     */
//    constexpr int r = 14;
//    /**
//     * Количество единиц в строках проверочной матрицы LDPC кода. Малая величина.
//     */
//    constexpr int w = 4;
//    /**
//     * Минимальное количество единиц в столбцах проверочной матрицы.
//     */
//    const int threshold_H = 1; // Не надо менять с единицы!
//    /**
//     * Минимальное количество единиц в столбцах систематической проверочной матрицы, кроме единичных столбцов.
//     * Не более r/2; на практике ~r/4.
//     */
//    const int threshold_H_sys = 3;
//    std::set<std::vector<int>> unique_cols; // Набор уникальных столбцов.
//    std::vector<int> column;
//    for (;;) {
//       H.clear();
//       for (int i = 0; i < r; ++i)
//          H.push_back(rng_n::get_random_row<n, w>(gd, ud));
//       bool is_ok = true;
//       unique_cols.clear();
//       for (int i = 0; i < n; ++i) {
//          int weight = 0;
//          column.clear();
//          for (int k = 0 ; k < r; ++k) {
//             weight += H.at(k).at(i);
//             column.push_back(H.at(k).at(i));
//          }
//          unique_cols.insert(column);
//          if (weight < threshold_H) {
//             is_ok = false;
//             break;
//          }
//       }
//       if (unique_cols.size() < (n - 2)) { // Подбирается, постепенно поднимая порог.
//          is_ok = false;
//       }
//       if (!is_ok)
//          continue;
//       std::tie(Hsys, swaps) = hamming::MakeParityMatrixSystematic(H, cyndrome, is_ok);
//       for (int i = 0; i < n - r; ++i) {
//          int weight = 0;
//          for (int k = 0; k < r; ++k)
//             weight += Hsys.at(k).at(i);
//          if (weight < threshold_H_sys) {
//             is_ok = false;
//             break;
//          }
//       }
//       if (is_ok) {
//          std::cout << "US: " << unique_cols.size() << std::endl;
//          break;
//       }
//    }
//    hamming::show_matrix(H, "H:");
//    hamming::show_matrix(Hsys, "Systematic form of H:");
//    hamming::show_vector(swaps, "Swaps:");
//    std::cout << "Threshold Hys: " << threshold_H_sys << ", r: " << r << ", w: " << w << std::endl;
// }


int main() {
   TestHamming();
   // MakeLdpc();
   return 0;
}