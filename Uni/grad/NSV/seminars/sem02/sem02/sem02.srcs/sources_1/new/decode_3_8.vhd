library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity decode_3_8 is
    Port ( a : in STD_LOGIC_VECTOR (2 downto 0);
           b : out STD_LOGIC_VECTOR (7 downto 0);
           en_n : in STD_LOGIC);
end decode_3_8;

architecture Behavioral of decode_3_8 is

begin

-- Dummy (but nicely readable) way
b <= "00000001" when a = "000" and en_n = '0' else
     "00000010" when a = "001" and en_n = '0' else
     "00000100" when a = "010" and en_n = '0' else
     "00001000" when a = "011" and en_n = '0' else
     "00010000" when a = "100" and en_n = '0' else
     "00100000" when a = "101" and en_n = '0' else
     "01000000" when a = "110" and en_n = '0' else
     "10000000" when a = "111" and en_n = '0' else
     "ZZZZZZZZ";

-- Using the process statement
--process(a) begin
--    case(a) is
--        when "000" => b <= "00000001";
--        when "001" => b <= "00000010";
--        when "010" => b <= "00000100";
--        when "011" => b <= "00001000";
--        when "100" => b <= "00010000";
--        when "101" => b <= "00100000";
--        when "110" => b <= "01000000";
--        when "111" => b <= "10000000";
--    end case;
--end process;

-- Usign the with-select statement
--with a select
--    b <= "00000001" when "000",
--         "00000010" when "001",
--         "00000100" when "010",
--         "00001000" when "011",
--         "00010000" when "100",
--         "00100000" when "101",
--         "01000000" when "110",
--         "10000000" when "111",
--         "00000000" when others;

end Behavioral;
