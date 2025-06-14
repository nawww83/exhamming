#include <iostream> // std::cout
#include "hamming.hpp"


static void TestHamming()
{
   using namespace hamming;
   using u8 = uint8_t;  // Тип внутреннего символа. В данном случае - байт.
   constexpr int R = 5; // Количество проверочных кодовых символов.
   constexpr int M = 4; // Количество внутренних символов на один кодовый символ (векторность).
   HammingExtended< u8, R, M > code;
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
   show_codeword(s, code.K, "Codeword:"); 
   // Вычисление синдрома.
   const auto& c = code.CalcSyndrome( s );
   show_cyndrome(c, "Cyndrome of the codeword:"); 
   // Стирание нескольких символов в канале. Кодовое расстояние кода равно 4,
   // кратность стирания q <= d - 1 = 3. Значит все стирания из не более трех 
   // символов будут успешно восстановлены. Но иногда восстанавливает и более высокие кратности.
   auto v = s;
   v[ 5 ] = { .mStatus = SymbolStatus::Erased, .mSymbol = {} };
   v[ 7 ] = { .mStatus = SymbolStatus::Erased, .mSymbol = {} };
   v[ 12 ] = { .mStatus = SymbolStatus::Erased, .mSymbol = {} };
   v[ 0 ] = { .mStatus = SymbolStatus::Erased, .mSymbol = {} };
   // Декодирование (восстановление стертых символов).
   int erased;
   [[maybe_unused]] int was_changed_strategy;
   const auto decode_is_ok = code.Decode( v, erased, was_changed_strategy );
   const bool recover_is_ok = v == a;
   show_codeword(v, code.K, "Decoded symbols:");  
   std::cout << "Decode is : " << (decode_is_ok ? "Ok\n" : "Failed\n");
   std::cout << "Recover is : " << (recover_is_ok ? "Ok\n" : "Failed\n");
   std::cout << '\n';
}


int main() {
   TestHamming();
   return 0;
}